import json
import os
import glob
import argparse
import logging
from typing import Dict, Any, Optional, List
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

DEFAULT_TEMPLATE_FOLDER = "templates"
DEFAULT_OUTPUT_FOLDER = "output"
LANGUAGES = ["EN", "CN"]

def load_data(file_path: str) -> Optional[Dict[str, Any]]:
    """Loads JSON data from the specified file path."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return data
    except (json.JSONDecodeError, FileNotFoundError, IOError) as e:
        logger.error(f"Failed to load data from {file_path}: {e}")
        return None

def determine_template_name(is_couple: bool, lang: str) -> str:
    """Determines the appropriate template filename based on couple status and language."""
    if is_couple:
        return f"carrier_couples_{lang}.html"
    return f"carrier_{lang}.html"

def generate_pdf(
    html_content: str, 
    output_path: str, 
    base_url: str = '.'
) -> bool:
    """Generates a PDF from HTML content using WeasyPrint."""
    try:
        HTML(string=html_content, base_url=base_url).write_pdf(output_path)
        return True
    except Exception as e:
        logger.error(f"Failed to generate PDF at {output_path}: {e}")
        return False

def process_file(
    file_path: str, 
    env: Environment, 
    output_folder: str, 
    languages: List[str],
    base_url: str = '.' # Added base_url parameter
) -> None:
    """Processes a single JSON file to generate reports in all specified languages."""
    data = load_data(file_path)
    if not data:
        return

    # Strict validation: data must be a dictionary (not a list like batch analysis)
    if not isinstance(data, dict):
        logger.debug(f"Skipping {file_path}: Content is not a dictionary (likely batch analysis).")
        return

    # Basic validations
    if 'report_metadata' not in data and 'primary_patient' not in data:
         logger.warning(f"File {file_path} does not appear to be a valid report data file. Skipping.")
         return

    logger.info(f"Processing file: {file_path}")

    is_couple = "partner" in data and data["partner"] is not None
    order_id = data.get('report_metadata', {}).get('order_id', 'Unknown')
    raw_name = data.get('primary_patient', {}).get('name', 'Patient')
    patient_name = "".join([c if c.isalnum() else "_" for c in raw_name]).replace('__', '_')

    for lang in languages:
        template_name = determine_template_name(is_couple, lang)
        
        try:
            template = env.get_template(template_name)
            html_out = template.render(data=data)
            
            file_name = f"Report_{order_id}_{patient_name}_{lang}.pdf"
            output_path = os.path.join(output_folder, file_name)
            
            if generate_pdf(html_out, output_path, base_url=base_url):
                logger.info(f"   ✅ Generated {lang}: {file_name}")
            else:
                 logger.error(f"   ❌ Failed to generate PDF for {lang}")

        except Exception as e:
            logger.error(f"   ❌ Error processing language {lang} for {file_path}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Generate PDF reports from JSON data.")
    parser.add_argument("--templates", default=DEFAULT_TEMPLATE_FOLDER, help="Path to templates folder")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_FOLDER, help="Path to output folder")
    parser.add_argument("--input", default=".", help="Path to input folder with JSON files")
    
    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.output):
        try:
            os.makedirs(args.output)
            logger.info(f"Created output directory: {args.output}")
        except OSError as e:
            logger.critical(f"Could not create output directory {args.output}: {e}")
            return

    # Check templates directory
    if not os.path.exists(args.templates):
        logger.critical(f"Templates directory not found: {args.templates}")
        return

    # Configure Jinja2
    env = Environment(loader=FileSystemLoader(args.templates))

    # Find JSON files
    input_pattern = os.path.join(args.input, "*.json")
    json_files = glob.glob(input_pattern)
    
    if not json_files:
        logger.warning(f"No JSON files found in {args.input}")
        return

    for json_file in json_files:
        process_file(json_file, env, args.output, LANGUAGES)

if __name__ == "__main__":
    main()
