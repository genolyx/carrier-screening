import os
import sqlite3
from datetime import datetime
from flask import Flask, render_template, request, redirect, url_for, jsonify
from google import genai
from google.genai import errors
from dotenv import load_dotenv
from pydantic import BaseModel
from flask_cors import CORS
import sys
import json

# Add Carrier Result path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Carrier_result')))
try:
    import generate as carrier_generate
except ImportError:
    print("Warning: Could not import Carrier_result.generate")
    carrier_generate = None

load_dotenv()
client = genai.Client(api_key=os.getenv("GEMINI_API_KEY"))

app = Flask(__name__)
CORS(app) # Enable CORS for all routes (needed for Portal to POST here)
DATABASE = 'genetic_knowledge.db'
MODEL_NAME = "gemini-2.5-flash-lite"

# Import Refactored Classifier
from acmg_classifier import ACMGClassifier
classifier = ACMGClassifier(model_name=MODEL_NAME)

class GeneKnowledge(BaseModel):
    omim_number: str
    function_summary: str
    function_summary: str
    disease_association: str
    function_summary: str
    disease_association: str
    inheritance: str
    disorder_name: str



# --- DATABASE SETUP ---
def init_db():
    """Initializes the database with expanded columns for Gene and Variant data."""
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        # Gene table now includes OMIM number and split knowledge sections
        cursor.execute('''CREATE TABLE IF NOT EXISTS gene_data (
                            gene_symbol TEXT PRIMARY KEY, 
                            function_summary TEXT, 
                            disease_association TEXT,
                            omim_number TEXT,
                            inheritance TEXT,
                            disorder TEXT)''')
        # Variant table now links to Gene table and includes ACMG classification
        cursor.execute('''CREATE TABLE IF NOT EXISTS variant_data (
                            id INTEGER PRIMARY KEY AUTOINCREMENT,
                            gene_symbol TEXT,
                            hgvs_result TEXT,
                            analysis_text TEXT,
                            acmg_classification TEXT,
                            variant_summary TEXT,
                            acmg_evidence TEXT,
                            FOREIGN KEY(gene_symbol) REFERENCES gene_data(gene_symbol) ON DELETE CASCADE)''')
        conn.commit()

init_db()

# --- HELPER FUNCTIONS ---
def get_gene_knowledge_with_db(mutation_input):
    """Retrieves Gene data from DB if available; otherwise, queries Gemini."""
    gene_symbol = mutation_input.split(' ')[0].upper().strip()
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT function_summary, disease_association, omim_number, inheritance, disorder FROM gene_data WHERE gene_symbol = ?", (gene_symbol,))
        row = cursor.fetchone()
        if row and row[0]:
            return f"OMIM: {row[2]}\nDISORDER: {row[4]}\nINHERITANCE: {row[3]}\nFUNCTION: {row[0]}\nDISEASE: {row[1]}"

    try:
        # STEP 1: Search for information (Tools enabled, returns Text)
        search_prompt = f"Using Google Search, find the latest technical details and OMIM data for the gene {gene_symbol}. Prioritize sources: OMIM, GeneReviews, UniProt, GTEx. Summarize the findings comprehensively. Also identify the specific name of the associated disorder (e.g. MTHFR deficiency)."
        
        search_response = client.models.generate_content(
            model=MODEL_NAME, 
            contents=search_prompt,
            config={'tools': [{'google_search': {}}]}
        )
        found_text = search_response.text

        # STEP 2: Extract structured data from the search results (No tools, returns JSON)
        extraction_prompt = f"""
        Extract the following details from the text below into the specified JSON format.
        TEXT: {found_text}
        
        Note for 'disorder_name': Provide the specific name of the primary disorder associated with this gene.
        """

        response = client.models.generate_content(
            model=MODEL_NAME, 
            contents=extraction_prompt,
            config={
                'response_mime_type': 'application/json',
                'response_schema': GeneKnowledge
            }
        )
        
        # Parse the structured response
        gene_data = response.parsed
        omim_part = gene_data.omim_number
        func_part = gene_data.function_summary
        disease_part = gene_data.disease_association
        inheritance_part = gene_data.inheritance
        disorder_part = gene_data.disorder_name

        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            cursor.execute("INSERT OR REPLACE INTO gene_data VALUES (?,?,?,?,?,?)", (gene_symbol, func_part, disease_part, omim_part, inheritance_part, disorder_part))
            conn.commit()
        return f"OMIM: {omim_part}\nDISORDER: {disorder_part}\nINHERITANCE: {inheritance_part}\nFUNCTION: {func_part}\nDISEASE: {disease_part}"
    except Exception as e:
        print(f"Error in get_gene_knowledge_with_db: {e}")
        return ""

def get_variant_analysis_with_db(mutation_input, context=None):
    """Retrieves Variant analysis from DB or Gemini to save quota."""
    parts = mutation_input.strip().split(' ', 1)
    if len(parts) < 2: return ""
    
    gene_symbol = parts[0].strip().upper()
    variant_hgvs = parts[1].strip()

    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT analysis_text, acmg_evidence FROM variant_data WHERE gene_symbol = ? AND hgvs_result = ?", (gene_symbol, variant_hgvs))
        row = cursor.fetchone()
        
        # If we have data, we might return it. 
        # BUT if user provides fresh context that might change classification (like AF), should we re-analyze?
        # For now, let's trust the cache to avoid quota burn, unless 'force_refresh' is implemented.
        # However, to be safe, if we have context, maybe we append it to the stored summary if missing?
        # Let's start by using it for NEW analysis only.
        if row and row[0]: return row[0]

    try:
        # Use Refactored Classifier
        data = classifier.classify(gene_symbol, variant_hgvs, context)
        
        text = data.analysis_text
        acmg = data.acmg_classification
        evidence = data.acmg_evidence

        # Separate Summary Generation (Creative Writing)
        summary_prompt = (f"Write a 3-4 sentence comprehensive molecular summary for the genetic variant {gene_symbol} {variant_hgvs}. "
                          f"Use the following analysis context: {text}. "
                          f"The consensus classification is {acmg}. "
                          f"CRITICAL: If gnomAD Allele Frequency is mentioned in the text, you MUST explicitly state it in this summary (e.g. 'Present in gnomAD at 0.04%').")
        
        summary_resp = client.models.generate_content(model=MODEL_NAME, contents=summary_prompt)
        summary = summary_resp.text

        # SAVE TO REPORTER DB
        # Convert List to JSON string for storage
        evidence_json = json.dumps(evidence)
        with sqlite3.connect(DATABASE) as conn:
            cursor = conn.cursor()
            cursor.execute("INSERT INTO variant_data (gene_symbol, hgvs_result, analysis_text, acmg_classification, variant_summary, acmg_evidence) VALUES (?,?,?,?,?,?)", 
                           (gene_symbol, variant_hgvs, text, acmg, summary, evidence_json))
            conn.commit()
            
        # [NEW] SYNC BACK TO PORTAL DB
        try:
             # Define Path relative to this file
             portal_db_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../phenotype_portal/variants.db'))
             if os.path.exists(portal_db_path):
                 with sqlite3.connect(portal_db_path) as pconn:
                     pcursor = pconn.cursor()
                     # Check if exists to update or insert
                     pcursor.execute("SELECT id FROM curated_variants WHERE gene = ? AND hgvsc = ?", (gene_symbol, variant_hgvs))
                     existing = pcursor.fetchone()
                     
                     source_note = "Auto-classified by Genetic Reporter AI"
                     if existing:
                         pcursor.execute("UPDATE curated_variants SET classification = ?, source = ?, notes = ? WHERE id = ?", 
                                        (acmg, "Genetic Reporter", source_note, existing[0]))
                     else:
                         pcursor.execute("INSERT INTO curated_variants (gene, hgvsc, hgvsp, classification, source, notes) VALUES (?, ?, ?, ?, ?, ?)",
                                        (gene_symbol, variant_hgvs, context.get('p') if context else "", acmg, "Genetic Reporter", source_note))
                     pconn.commit()
                 print(f"Successfully synced classification '{acmg}' to Portal DB.")
        except Exception as e:
            print(f"Failed to sync with Portal DB: {e}")

        return text
    except Exception as e:
        print(f"Error in get_variant_analysis_with_db: {e}")
        return ""

# --- ROUTES ---
@app.route('/')
def index():
    # Pass dummy patient data so the template doesn't crash on load
    dummy_patient = {'name': 'Guest', 'age': 'N/A', 'sex': 'N/A', 'notes': 'No patient loaded.'}
    return render_template('report.html', p=dummy_patient, variants=[])

@app.route('/test_acmg', methods=['GET', 'POST'])
def test_acmg():
    result = None
    if request.method == 'POST':
        gene = request.form.get('gene')
        hgvs = request.form.get('hgvs')
        context = {
            'af': request.form.get('af'),
            'clinvar': request.form.get('clinvar'),
            'effect': request.form.get('effect'),
            'chrom': request.form.get('chrom'),
            'pos': request.form.get('pos'),
            'ref': request.form.get('ref'),
            'alt': request.form.get('alt'),
            # Default Stubs
            'nm': 'NM_TEST',  
            'p': 'p.Test',
            'rsid': 'rs000'
        }
        
        try:
            # Direct call to classifier logic
            result = classifier.classify(gene, hgvs, context)
        except Exception as e:
            result = {"acmg_classification": "ERROR", "acmg_evidence": str(e), "analysis_text": ""}

    return render_template('test_classifier.html', result=result)

@app.route('/genes')
def list_genes():
    """Displays the split-view Knowledge Library Manager."""
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM gene_data ORDER BY gene_symbol ASC")
        genes = cursor.fetchall()
        cursor.execute("SELECT * FROM variant_data ORDER BY gene_symbol ASC, hgvs_result ASC")
        variants = cursor.fetchall()
    return render_template('genes.html', genes=genes, variants=variants)

@app.route('/delete_gene/<symbol>')
def delete_gene(symbol):
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("DELETE FROM gene_data WHERE gene_symbol = ?", (symbol,))
        conn.commit()
    return redirect(url_for('list_genes'))

@app.route('/delete_variant/<int:vid>')
def delete_variant(vid):
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        cursor.execute("DELETE FROM variant_data WHERE id = ?", (vid,))
        conn.commit()
    return redirect(url_for('list_genes'))

@app.route('/update_variant', methods=['POST'])
def update_variant():
    """Updates manually edited variant data AND gene data."""
    vid = request.form.get('vid')
    classification = request.form.get('classification')
    summary = request.form.get('summary')
    evidence = request.form.get('evidence')
    
    # New Gene Fields
    omim = request.form.get('omim')
    inheritance = request.form.get('inheritance')
    disorder = request.form.get('disorder')
    disease_assoc = request.form.get('disease_association')

    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        
        # 1. Update Variant Data
        cursor.execute("""
            UPDATE variant_data 
            SET acmg_classification = ?, variant_summary = ?, acmg_evidence = ?
            WHERE id = ?
        """, (classification, summary, evidence, vid))
        
        # 2. Update Gene Data (via gene_symbol lookup)
        # First get the gene symbol for this variant
        cursor.execute("SELECT gene_symbol FROM variant_data WHERE id = ?", (vid,))
        row = cursor.fetchone()
        if row:
            gene_symbol = row[0]
            cursor.execute("""
                UPDATE gene_data 
                SET omim_number = ?, inheritance = ?, disorder = ?, disease_association = ?
                WHERE gene_symbol = ?
            """, (omim, inheritance, disorder, disease_assoc, gene_symbol))
            
        conn.commit()
    
    if request.args.get('ajax'):
        return jsonify({"status": "success", "message": "Variant updated successfully"})
        
    return redirect(url_for('list_genes'))

@app.route('/export_json')
def export_json():
    """Exports gene and variant data to a JSON file."""
    with sqlite3.connect(DATABASE) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        
        # Query joining Gene and Variant tables
        query = """
            SELECT 
                g.gene_symbol as gene,
                g.omim_number,
                g.inheritance,
                g.disorder,
                v.hgvs_result as mutation,
                v.acmg_classification as classification,
                v.acmg_evidence,
                g.disease_association as gene_description,
                g.function_summary as gene_function,
                v.variant_summary as variant_summary,
                v.analysis_text as full_analysis
            FROM variant_data v
            JOIN gene_data g ON v.gene_symbol = g.gene_symbol
            ORDER BY g.gene_symbol, v.hgvs_result
        """
        cursor.execute(query)
        rows = cursor.fetchall()
        
        # Convert rows to list of dicts and format fields
        data = []
        for row in rows:
            item = dict(row)
            # Enforce Title Case for specific fields as requested
            if item.get('inheritance'):
                item['inheritance'] = item['inheritance'].title()
            if item.get('classification'):
                item['classification'] = item['classification'].title()
            data.append(item)
        
    # Create response with JSON file attachment
    from flask import Response
    import json
    return Response(
        json.dumps(data, indent=2),
        mimetype="application/json",
        headers={"Content-Disposition": "attachment;filename=genetic_library_export.json"}
    )

@app.route('/batch_analyze', methods=['GET', 'POST'])
def batch_analyze():
    """Handles batch input, triggers analysis, and returns a JSON export."""
    if request.method == 'GET':
        return render_template('batch.html')
    
    # POST request
    raw_input = request.form.get('batch_input', '')
    lines = [line.strip() for line in raw_input.split('\n') if line.strip()]
    
    processed_results = []
    
    for line in lines:
        try:
            # 1. Parse Input
            # Expecting "GENE VARIANT" e.g. "MTHFR c.677C>T"
            parts = line.split(' ', 1)
            if len(parts) < 2:
                continue # Skip malformed lines
                
            gene_symbol_in = parts[0].strip().upper()
            variant_hgvs_in = parts[1].strip()
            full_mutation_str = f"{gene_symbol_in} {variant_hgvs_in}"
            
            # 2. Trigger Analysis (this functions update the DB)
            # Fetch Gene Data (Step 1 of cache-aside logic)
            get_gene_knowledge_with_db(full_mutation_str)
            # Fetch Variant Data (Step 2 of cache-aside logic)
            get_variant_analysis_with_db(full_mutation_str)
            
            # 3. Retrieve from DB for Export
            with sqlite3.connect(DATABASE) as conn:
                conn.row_factory = sqlite3.Row
                cursor = conn.cursor()
                query = """
                    SELECT 
                        g.gene_symbol as gene,
                        g.omim_number,
                        g.inheritance,
                        g.disorder,
                        v.hgvs_result as mutation,
                        v.acmg_classification as classification,
                        v.acmg_evidence,
                        g.disease_association as gene_description,
                        g.function_summary as gene_function,
                        v.variant_summary as variant_summary,
                        v.analysis_text as full_analysis
                    FROM variant_data v
                    JOIN gene_data g ON v.gene_symbol = g.gene_symbol
                    WHERE g.gene_symbol = ? AND v.hgvs_result = ?
                """
                cursor.execute(query, (gene_symbol_in, variant_hgvs_in))
                row = cursor.fetchone()
                if row:
                    item = dict(row)
                    # Enforce Title Case for specific fields as requested
                    if item.get('inheritance'):
                        item['inheritance'] = item['inheritance'].title()
                    if item.get('classification'):
                        item['classification'] = item['classification'].title()
                    processed_results.append(item)
                    
        except Exception as e:
            print(f"Error processing batch item {line}: {e}")
            continue

    # Return JSON file
    from flask import Response
    import json
    
    # Check for clinical synthesis (from report page export)
    synthesis_text = request.form.get('clinical_synthesis')
    
    if synthesis_text:
        final_output = {
            "clinical_interpretation": synthesis_text,
            "variants": processed_results
        }
    else:
        # Maintain backward compatibility for raw batch analysis
        final_output = processed_results

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Response(
        json.dumps(final_output, indent=2),
        mimetype="application/json",
        headers={"Content-Disposition": f"attachment;filename=batch_analysis_{timestamp}.json"}
    )

@app.route('/generate_batch', methods=['GET', 'POST'])
def generate_batch():
    if request.method == 'GET':
        return "<h1>Genetic Reporter API</h1><p>This endpoint receives data from the Phenotype Portal.</p><p>Please use the Portal to select variants and click 'Send to Reporter'.</p>", 200
    """Generates a consolidated clinical report for MULTIPLE variants."""
    patient_data = {
        "name": request.form.get('patient_name'),
        "age": request.form.get('age'),
        "sex": request.form.get('sex'),
        "notes": request.form.get('clinical_notes')
    }
    
    raw_input = request.form.get('batch_input', '')
    
    # Try parsing as JSON first (new rich format)
    import json
    parsed_variants = []
    try:
        parsed_variants = json.loads(raw_input)
        if not isinstance(parsed_variants, list):
            parsed_variants = []
    except:
        pass

    if not parsed_variants and raw_input:
        # Fallback to newline splitting (legacy string format)
        lines = [line.strip() for line in raw_input.split('\n') if line.strip()]
        for line in lines:
            parts = line.split(' ', 1)
            if len(parts) >= 2:
                parsed_variants.append({
                    'gene': parts[0].strip().upper(),
                    'hgvsc': parts[1].strip(),
                    'original_text': line
                })

    if not parsed_variants:
        return "No variants provided", 400

    combined_background = ""
    findings_list = []
    variants_data = []

    for v_obj in parsed_variants:
        # Extract enriched context
        gene_sym = v_obj.get('gene')
        hgvsc = v_obj.get('hgvsc')
        full_mutation_str = v_obj.get('original_text', f"{gene_sym} {hgvsc}")
        
        # Build context dict for analysis
        context = {
            'nm': v_obj.get('transcript'),
            'p': v_obj.get('protein'),
            'effect': v_obj.get('effect'),
            'af': v_obj.get('af'),
            'rsid': v_obj.get('rsid'),
            'clinvar': v_obj.get('clinvar'),
            'portal_class': v_obj.get('acmg_class'),
            'chrom': v_obj.get('chrom'),
            'pos': v_obj.get('pos'),
            'ref': v_obj.get('ref'),
            'alt': v_obj.get('alt')
        }
        
        findings_list.append(full_mutation_str)

        # Trigger Analysis (Ensure DB is populated) with enriched context
        gene_info = get_gene_knowledge_with_db(full_mutation_str)
        var_text_analysis = get_variant_analysis_with_db(full_mutation_str, context=context)
        
        # Retrieve structured data for template
        with sqlite3.connect(DATABASE) as conn:
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            c.execute("""
                SELECT 
                    v.id, 
                    v.variant_summary, 
                    v.acmg_evidence, 
                    v.acmg_classification,
                    g.disease_association,
                    g.function_summary,
                    g.omim_number,
                    g.inheritance,
                    g.disorder
                FROM variant_data v
                LEFT JOIN gene_data g ON v.gene_symbol = g.gene_symbol
                WHERE v.gene_symbol = ? AND v.hgvs_result = ?
            """, (gene_sym, hgvsc))
            row = c.fetchone()

        if row:
            # Pre-process evidence to ensure line breaks between codes
            # Looks for patterns like "PM2:", "PVS1:", "PP3:" preceded by whitespace/punctuation
            raw_evidence = row['acmg_evidence']
            formatted_evidence = "Evidence pending..."
            
            # Try parsing as JSON List (New 2.0 format)
            try:
                import json
                parsed_ev = json.loads(raw_evidence)
                if isinstance(parsed_ev, list):
                    # Fix: Stitch broken lines (legacy data support)
                    stitched_ev = []
                    for ev in parsed_ev:
                        clean = str(ev).replace('\n', ' ').strip()
                        if not clean: continue
                        
                        is_new_item = False
                        upper = clean.upper()
                        # Strict Prefix Check
                        if upper.startswith("CLINVAR") or upper.startswith("NOTE"):
                            is_new_item = True
                        elif len(clean) > 3 and clean[0:2].isalpha() and clean[2].isdigit() and clean[3] == ':':
                             is_new_item = True
                        elif clean.startswith("BA1") or clean.startswith("BS1"):
                             is_new_item = True
                        
                        if is_new_item or not stitched_ev:
                            stitched_ev.append(clean)
                        else:
                            stitched_ev[-1] += " " + clean # Merge with previous

                    formatted_evidence = "\n".join(stitched_ev)
                else:
                    raise ValueError("Not a list")
            except:
                # Fallback to Legacy Regex/String handling
                if raw_evidence:
                    import re
                    # Replace semicolons or comma-space sequences with newline
                    formatted_evidence = re.sub(r'[;,]\s*', '\n', str(raw_evidence))
                    # Ensure separate lines for codes that might just be space-separated
                    formatted_evidence = re.sub(r'\s+(?=[A-Z]{2,4}\d+:)', '\n', formatted_evidence)

            variants_data.append({
                'gene': gene_sym,
                'mutation': hgvsc,
                'summary': row['variant_summary'] or "Summary pending...",
                'evidence': formatted_evidence,
                'classification': row['acmg_classification'] or "Uncertain Significance",
                'gene_description': row['disease_association'] or "Description pending...",
                'gene_function': row['function_summary'] or "",
                'omim': row['omim_number'] or "",
                'inheritance': row['inheritance'] or "",
                'disorder': row['disorder'] or "",
                'id': row['id'],
                'gene_info': gene_info
            })
            
            # Context for Synthesis Prompt
            combined_background += f"\n--- VARIANT: {full_mutation_str} ---\n"
            combined_background += f"CLASSIFICATION: {row['acmg_classification']}\n"
            combined_background += f"SUMMARY: {row['variant_summary']}\n"
            combined_background += f"GENE INFO: {gene_info}\n"

    findings_str = ", ".join(findings_list)

    # Simplified Prompt: ONLY Clinical Interpretation
    # We will render the specific variant details using Jinja from the DB data directly
    report_prompt = f"""
    You are a Senior Clinical Geneticist. Write a 'Clinical Interpretation' section for a patient report.
    
    PATIENT: {patient_data['name']}
    AGE/SEX: {patient_data['age']} / {patient_data['sex']}
    FINDINGS: {findings_str}
    CLINICAL CONTEXT: {patient_data['notes']}
    
    The detailed evidence and molecular summaries are already prepared. 
    YOUR JOB: Write a cohesive "Clinical Interpretation" paragraph (or two) that synthesizes these findings.
    - Do they explain the user's symptoms?
    - Are there synergistic effects?
    - What is the overall clinical impression?
    
    CONTEXT DATA:
    {combined_background}
    
    FORMATTING: Use HTML <p> tags. Do NOT use headers like 'Molecular Findings' or 'ACMG', those are already handled. Just the synthesis text.
    """

    try:
        response = client.models.generate_content(model=MODEL_NAME, contents=report_prompt)
        synthesis_text = response.text
    except Exception as e:
        print(f"Error in generate_batch report: {e}")
        synthesis_text = "<p>API Error generating synthesis.</p>"

    # Render Report using STRUCTRED DATA
    patient_data['mutation'] = findings_str 
    patient_data['raw_input'] = raw_input # Pass original input for JSON export button
    final_html = render_template('report.html', 
                                synthesis=synthesis_text, 
                                variants=variants_data, 
                                p=patient_data)

    # Save Report
    if not os.path.exists('reports'):
        os.makedirs('reports')
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
    safe_name = "".join([c if c.isalnum() else "_" for c in patient_data['name']]).lower()
    filename = f"reports/batch_report_{safe_name}_{timestamp}.html"
    with open(filename, "w") as f:
        f.write(final_html)
    
    return final_html

@app.route('/download_report/<path:filename>')
def download_report(filename):
    """Serves the generated PDF report."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    carrier_dir = os.path.join(base_dir, "..", "Carrier_result", "output")
    from flask import send_from_directory
    return send_from_directory(carrier_dir, filename)

@app.route('/generate_pdf_carrier', methods=['POST'])
def generate_pdf_carrier():
    if not carrier_generate:
        return jsonify({"status": "error", "message": "Carrier generation module not loaded"}), 500
        
    try:
        data = request.json
        if not data:
             return jsonify({'status': 'error', 'message': 'No JSON data provided'}), 400

        # Setup Paths
        base_dir = os.path.dirname(os.path.abspath(__file__))
        carrier_dir = os.path.join(base_dir, "..", "Carrier_result")
        template_dir = os.path.join(carrier_dir, "templates")
        output_dir = os.path.join(carrier_dir, "output")
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Standard Jinja Env
        from jinja2 import Environment, FileSystemLoader
        env = Environment(loader=FileSystemLoader(template_dir))
        
        # Save Request Data to Persistent Folder
        json_output_dir = os.path.join(carrier_dir, "input") # Use 'input' as that's standard for generate.py
        if not os.path.exists(json_output_dir):
            os.makedirs(json_output_dir)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        safe_name = "".join([c if c.isalnum() else "_" for c in data.get('primary_patient', {}).get('name', 'Patient')]).replace('__', '_')
        json_filename = f"Input_{safe_name}_{timestamp}.json"
        json_path = os.path.join(json_output_dir, json_filename)
        
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)

        # Generate BOTH languages using the shared module
        # process_file(file_path, env, output_folder, languages, base_url)
        carrier_generate.process_file(json_path, env, output_dir, ["EN", "CN"], base_url=carrier_dir)
        
        # Kept the JSON file (no cleanup)
        # if os.path.exists(temp_path):
        #    os.remove(temp_path)
            
        return jsonify({
            'status': 'success', 
            'message': f'Generated EN & CN reports in {output_dir}'
        })
        
    except Exception as e:
        print(f"Error generating PDF: {e}")
        return jsonify({'status': 'error', 'message': str(e)}), 500

@app.route('/manual_add', methods=['POST'])
def manual_add():
    """Handles manual data entry to bypass API quotas."""
    entry_type = request.form.get('type')
    symbol = request.form.get('symbol').strip().upper()
    with sqlite3.connect(DATABASE) as conn:
        cursor = conn.cursor()
        if entry_type == 'gene':
            omim = request.form.get('omim').strip()
            func = request.form.get('function').strip()
            disease = request.form.get('disease').strip()
            # Default inheritance/disorder to empty for manual entry until UI is updated
            inheritance = ""
            disorder = ""
            cursor.execute("INSERT OR REPLACE INTO gene_data VALUES (?,?,?,?,?,?)", (symbol, func, disease, omim, inheritance, disorder))
        else:
            content = request.form.get('content').strip()
            # Split input symbol into Gene and Variant manually if needed or assume symbol is Gene and add separate field?
            # For "Bypass AI Quota", user likely enters "Gene Variant".
            parts = symbol.split(' ', 1)
            if len(parts) == 2:
                g_sym = parts[0].strip().upper()
                v_hgvs = parts[1].strip()
                cursor.execute("INSERT INTO variant_data (gene_symbol, hgvs_result, analysis_text) VALUES (?,?,?)", (g_sym, v_hgvs, content))
        conn.commit()
    return redirect(url_for('list_genes'))

@app.route('/generate', methods=['POST'])
def generate():
    """Generates the final clinical report using the preserved 5-section prompt."""
    data = {
        "name": request.form.get('patient_name'),
        "age": request.form.get('age'),
        "sex": request.form.get('sex'),
        "mutation": request.form.get('mutation'),
        "notes": request.form.get('clinical_notes')
    }

    gene_info = get_gene_knowledge_with_db(data['mutation'])
    var_info = get_variant_analysis_with_db(data['mutation'])
    combined = (gene_info or "Pending Gene Research") + "\n\n" + (var_info or "Pending Variant Analysis")

    # This prompt is preserved exactly as requested
    report_prompt = f"""
    You are a Senior Clinical Geneticist. Draft a formal clinical summary.
    
    PATIENT: {data['name']} | FINDING: {data['mutation']} | CONTEXT: {data['notes']}
    BACKGROUND DATA: {combined}
    
    FORMATTING: Use HTML <h3>, <p>, <strong>, <ul>, and <li> tags.
    
    CONTENT STRUCTURE:
    1. <h3>Molecular Finding:</h3> Describe the mutation and gene function.
    2. <h3>Clinical Interpretation:</h3> Link finding to patient symptoms ({data['notes']}).
    3. <h3>Gene & Disease Overview:</h3> General overview of disease and inheritance using BACKGROUND DATA.
    4. <h3>ACMG Interpretation Summary:</h3> 
       - Provide a bulleted list (<ul> and <li>) summarizing evidence (e.g., PM2, PP3) 
         specific to the variant {data['mutation']}.
    5. <h3>Classification:</h3> 
       - Based on the evidence above, provide the final classification (e.g., PATHOGENIC, LIKELY PATHOGENIC, or VARIANT OF UNCERTAIN SIGNIFICANCE). 
       - Format this as: <p><strong>Final Classification:</strong> [CLASSIFICATION HERE]</p>
    """

    try:
        response = client.models.generate_content(model=MODEL_NAME, contents=report_prompt)
        final_summary = response.text
    except Exception as e:
        print(f"Error in generate report: {e}")
        final_summary = f"<p>API Error. Displaying raw data:</p><p>{combined}</p>"

    final_html = render_template('report.html', summary=final_summary, p=data)

    if not os.path.exists('reports'):
        os.makedirs('reports')

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
    safe_name = "".join([c if c.isalnum() else "_" for c in data['name']]).lower()
    filename = f"reports/report_{safe_name}_{timestamp}.html"
    with open(filename, "w") as f:
        f.write(final_html)
    
    return final_html
@app.route('/batch_from_portal', methods=['POST'])
def batch_from_portal():
    """Receives selected variants from Phenotype Portal and pre-populates the DB."""
    data = request.json
    variants = data.get('variants', [])
    
    if not variants:
        return jsonify({"status": "error", "message": "No variants provided"}), 400

    count = 0
    for full_mutation_str in variants:
        try:
            # Trigger Analysis (populates the DB)
            get_gene_knowledge_with_db(full_mutation_str)
            get_variant_analysis_with_db(full_mutation_str)
            count += 1
        except Exception as e:
            print(f"Error analyzing {full_mutation_str}: {e}")

    return jsonify({"status": "success", "count": count, "message": f"Processed {count} variants"})


if __name__ == '__main__':
    app.run(debug=True, port=5001, host='0.0.0.0')
