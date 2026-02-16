"""
Script to scrape ATcT thermochemical data.

Usage:
    python scripts/scrape_atct.py -o atct_raw.csv

This will need to be customized based on the actual ATcT page structure.
You may need to use Selenium or BeautifulSoup to navigate the interactive page.
"""

import argparse
import pandas as pd

# TODO: Implement web scraping
# The ATcT page likely requires JavaScript interaction
# Suggested approach:
# 1. Use Selenium to navigate the page
# 2. Extract species names and ΔHf° values (298.15 K, shown in red)
# 3. Save to CSV with columns: name, formula, dhf_298K_kcal_mol, uncertainty

def scrape_atct():
    """Scrape ATcT database."""
    raise NotImplementedError(
        "Web scraping implementation needed. "
        "ATcT page requires JavaScript interaction. "
        "Consider using Selenium with ChromeDriver."
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default="atct_raw.csv")
    args = parser.parse_args()

    scrape_atct()
