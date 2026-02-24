from pydantic import BaseModel
from google import genai
import os
from dotenv import load_dotenv

load_dotenv()

# Pydantic Model for Structured Output
class VariantAnalysis(BaseModel):
    acmg_classification: str
    analysis_text: str
    acmg_evidence: list[str]

class ACMGClassifier:
    def __init__(self, model_name="gemini-2.5-flash-lite", api_key=None):
        self.api_key = api_key or os.getenv("GEMINI_API_KEY")
        # In Lite Mode (used by Portal), we might not have a key, and we don't need one.
        # Only crash if we try to use AI.
        if self.api_key:
            self.client = genai.Client(api_key=self.api_key)
        else:
            self.client = None
        self.model_name = model_name

    def _lookup_local_clinvar(self, chrom, pos, ref, alt):
        try:
            import pysam
            # Path relative to this file: ../phenotype_portal/clinvar/clinvar.vcf.gz
            vcf_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../phenotype_portal/clinvar/clinvar.vcf.gz'))
            if not os.path.exists(vcf_path): 
                return None
            
            # Simple chrom handling: try with and without 'chr'
            chroms_to_try = [chrom]
            if chrom.startswith('chr'): chroms_to_try.append(chrom[3:])
            else: chroms_to_try.append(f"chr{chrom}")
            
            with pysam.VariantFile(vcf_path) as vcf:
                for c in chroms_to_try:
                    try:
                        for rec in vcf.fetch(c, int(pos)-1, int(pos)+1):
                            if rec.pos == int(pos) and rec.ref == ref and alt in rec.alts:
                                # Extract CLNSIG
                                sig = rec.info.get('CLNSIG')
                                if sig:
                                    return sig[0] # Return first match
                    except ValueError:
                        continue # contig not found in VCF
        except Exception as e:
            print(f"Local ClinVar lookup warning: {e}")
        return None

    def classify(self, gene_symbol, variant_hgvs, context=None, lite_mode=False, do_local_lookup=True):
        """
        Orchestrates the ACMG classification process:
        1. Pre-flight deterministic checks (Frequency, ClinVar, Variant Type).
        2. (Optional) Local ClinVar Lookup (if do_local_lookup=True).
        3. AI Synthesis and Extraction (Skipped if lite_mode=True).
        """
        
        # 1. Build Context and Pre-flight Rules
        context_str = ""
        missing_data_instructions = ""
        known_facts = []
        
        if context:
            user_af = context.get('af')
            user_transcript = context.get('nm')
            user_clinvar = context.get('clinvar')
            user_effect = context.get('effect')
            user_portal_class = context.get('portal_class')
            
            # --- PRE-FLIGHT RULE LOGIC ---
            
            # 1. Frequency Rules
            if user_af is not None and str(user_af) not in ['None', '-', '']:
                try:
                    val = float(user_af)
                    if val > 0.05:
                        known_facts.append(f"FREQUENCY RULE: AF is {val:.4f} (>5%). This meets criterion BA1 (Stand-alone Benign).")
                        if lite_mode: return VariantAnalysis(acmg_classification="Benign", acmg_evidence=["BA1: Allele Frequency > 5%"], analysis_text="Pre-flight classification.")
                    elif val > 0.01:
                        known_facts.append(f"FREQUENCY RULE: AF is {val:.4f} (>1%). This meets criterion BS1 (Strong Benign).")
                        if lite_mode: return VariantAnalysis(acmg_classification="Likely Benign", acmg_evidence=["BS1: Allele Frequency > 1%"], analysis_text="Pre-flight classification.")
                    elif val < 0.0001:
                        # Note: 0.0 falls here, which is correct (Rare).
                        known_facts.append(f"FREQUENCY RULE: AF is {val:.6f} (Very Rare). This supports PM2 (Moderate Pathogenic).")
                except:
                    pass
            
            # 2. ClinVar Rules
            # Fallback: Local Lookup if missing
            if do_local_lookup and (not user_clinvar or user_clinvar in ['None', '-']) and context.get('chrom') and context.get('pos'):
                 found_sig = self._lookup_local_clinvar(context.get('chrom'), context.get('pos'), context.get('ref'), context.get('alt'))
                 if found_sig:
                     user_clinvar = found_sig
                     known_facts.append(f"LOCAL CLINVAR LOOKUP: {user_clinvar} (Found in local database).")

            if user_clinvar and user_clinvar not in ['None', '-']:
                known_facts.append(f"CLINVAR STATUS: {user_clinvar}. You must prioritize this classification.")
                # Case-insensitive checks
                cv_lower = user_clinvar.lower()
                if lite_mode:
                    if "pathogenic" in cv_lower:
                        # Catch "Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"
                        cls = "Pathogenic"
                        if "likely" in cv_lower and "pathogenic/likely" not in cv_lower:
                             cls = "Likely Pathogenic"
                        return VariantAnalysis(acmg_classification=cls, acmg_evidence=[f"ClinVar: {user_clinvar}"], analysis_text="Pre-flight classification.")
                    
                    if "benign" in cv_lower:
                        cls = "Benign"
                        if "likely" in cv_lower and "benign/likely" not in cv_lower:
                             cls = "Likely Benign"
                        return VariantAnalysis(acmg_classification=cls, acmg_evidence=[f"ClinVar: {user_clinvar}"], analysis_text="Pre-flight classification.")

            # 3. Effect Rules (PVS1 Candidate)
            eff_lower = str(user_effect).lower()
            if 'stop' in eff_lower or 'frameshift' in eff_lower or 'nonsense' in eff_lower:
                 known_facts.append(f"VARIANT TYPE: {user_effect}. This is a LoF variant. Check if gene is Haploinsufficient (PVS1).")
                 # We don't auto-classify PVS1 in lite mode without knowing gene constraint, unless we add gene constraint logic later.
                 # For now, return VUS or Potential.
        
        if lite_mode:
            # If we reached here in lite mode, we didn't hit a definitive Stop/Go rule.
            return VariantAnalysis(acmg_classification="VUS", acmg_evidence=["Lite Mode: No definitive rule hit."], analysis_text="Preliminary check only.")

            known_facts_str = "\n".join(known_facts)
            
            context_str = (f"CONTEXT PROVIDED BY USER:\n"
                           f"- Transcript: {user_transcript}\n"
                           f"- Protein: {context.get('p')}\n"
                           f"- Effect: {user_effect}\n"
                           f"- AF: {user_af}\n"
                           f"- dbSNP: {context.get('rsid')}\n"
                           f"- ClinVar: {user_clinvar}\n\n"
                           f"KNOWN FACTS (DO NOT HALLUCINATE):\n{known_facts_str}")
            

            if not user_transcript or user_transcript == 'None':
                 missing_data_instructions += " User did NOT provide specific transcript. Identify the canonical transcript (NM_...) using RefSeq/MANE."

        else:
            missing_data_instructions = " No user context provided. YOU MUST independenty search for Allele Frequency (gnomAD), Canonical Transcript, and Computational Predictions."

        # 2. Construct Search Prompt
        search_prompt = (f"Using Google Search, analyze the genetic variant {gene_symbol} {variant_hgvs}. "
                  f"{context_str} "
                  f"{missing_data_instructions} "
                  f"Find recent clinical interpretations, molecular impact, and ACMG criteria. "
                  f"Prioritize sources: ClinVar, ClinVarMiner, DECIPHER, ClinGen, HGMD, LOVD, gnomAD, DGV, "
                  f"PubMed, Mastermind, Varsome, Franklin (Genoox), MitoMap. "
                  f"MANDATORY LITERATURE SEARCH: explicitly search for '{gene_symbol} {variant_hgvs}' AND ('functional study' OR 'case report') on PubMed. "
                  f"If you find ANY papers describing this specific variant, you MUST cite the PMID in the summary or evidence. "
                  f"Specific Data Needed: 1. Computational predictions (REVEL, CADD, SpliceAI, PolyPhen-2, PROVEAN) - are they deleterious? "
                  f"3. Functional Evidence (PS3): Are there assays proving impact? Cite PMID. "
                  f"Explicitly find the consensus ACMG classification from ClinVar."
                  f"If ClinVar is inconclusive, use the found frequency and computational data to support codes like PM2 or PP3.")
        
        # 3. Call AI for Search
        search_response = self.client.models.generate_content(
            model=self.model_name, 
            contents=search_prompt,
            config={'tools': [{'google_search': {}}]}
        )
        found_text = search_response.text

        # 4. Construct Extraction Prompt
        extraction_prompt = f"""
        Extract details into the JSON format.
        TEXT: {found_text}
        
        INSTRUCTIONS FOR CLASSIFICATION:
        1. **Evaluate Criteria Explicitly**:
           - **PVS1**: Is it a Null variant (stop/frame) in a gene where LoF is a known mechanism? (Cite evidence).
           - **PM2/BA1**: Check gnomAD frequency. (Rare = PM2, Common >5% = BA1).
           - **PP3**: Check computational scores (REVEL > 0.75, CADD > 25).
           - **PS3**: Are there functional studies? (Cite PMID).
           - **ClinVar**: What is the consensus?
        
        2. **Synthesize Verdict**:
           - Use the ACMG Combination Rules (e.g. 1 Very Strong + 1 Strong = Pathogenic).
           - If ClinVar is "Pathogenic" and well-supported, align with it.
           
        Note for 'acmg_classification': The final calculated consensus (e.g. "Pathogenic", "VUS").
        Note for 'acmg_evidence': List the APPLIED codes with brief justification. Format: "Code: Description" (one per line).
        """

        # 5. Call AI for Extraction
        result = self.client.models.generate_content(
            model=self.model_name, 
            contents=extraction_prompt,
            config={
                'response_mime_type': 'application/json',
                'response_schema': VariantAnalysis
            }
        )
        
        parsed = result.parsed
        
        # 6. Post-Processing: Stitch broken evidence lines
        # AI sometimes returns ["PM2: Text...", "continued text"]
        final_evidence = []
        if parsed and parsed.acmg_evidence:
            for ev in parsed.acmg_evidence:
                clean = ev.replace('\n', ' ').replace('\r', '').strip()
                if not clean: continue
                
                # Check for standard prefixes: PVS1, PM2, PP3, BA1, BS1, ClinVar or standard code pattern
                is_new_item = False
                upper = clean.upper()
                
                # Known prefixes
                if upper.startswith("CLINVAR") or upper.startswith("NOTE"):
                    is_new_item = True
                # ACMG Code Pattern (e.g. PM2:, PVS1:)
                elif len(clean) > 3 and clean[0:2].isalpha() and clean[2].isdigit():
                     is_new_item = True
                # Catch "BA1" / "BS1"
                elif clean.startswith("BA1") or clean.startswith("BS1"):
                     is_new_item = True
                     
                if is_new_item:
                    final_evidence.append(clean)
                else:
                    # It's a continuation line (or garbage), append to previous if exists
                    if final_evidence:
                        final_evidence[-1] += " " + clean
                    else:
                        # Fallback if first line is weird
                        final_evidence.append(clean)
                        
            parsed.acmg_evidence = final_evidence
            
        return parsed
