from app2 import get_variant_analysis_with_db, init_db
import sqlite3

DATABASE = 'genetic_knowledge.db'

def test_pigq():
    print("Testing PIGQ c.560G>A analysis...")
    
    # 1. Clear existing data to force fresh Gemini call
    with sqlite3.connect(DATABASE) as conn:
        c = conn.cursor()
        c.execute("DELETE FROM variant_data WHERE gene_symbol='PIGQ' AND hgvs_result='c.560G>A'")
        conn.commit()
        
    # 2. Run analysis
    # Note: Requires correct GEMINI_API_KEY in .env
    result = get_variant_analysis_with_db("PIGQ c.560G>A")
    
    # 3. Check DB content
    with sqlite3.connect(DATABASE) as conn:
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        c.execute("SELECT * FROM variant_data WHERE gene_symbol='PIGQ' AND hgvs_result='c.560G>A'")
        row = c.fetchone()
        
        if row:
            print("\n--- RESULTS ---")
            print(f"Classification: {row['acmg_classification']}")
            print(f"Evidence:\n{row['acmg_evidence']}")
            print(f"Summary:\n{row['variant_summary']}")
        else:
            print("No data found.")

if __name__ == "__main__":
    test_pigq()
