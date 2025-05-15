from Bio import Entrez, SeqIO
import time
import os
import csv
import pandas as pd
import matplotlib.pyplot as pyplot
##Api key: 113dfe4392781f84489d30d39a15522f7208
class NCBIRetriever:
    def __init__(self, email, api_key,max_size,min_size):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key
        """Tworzenie zmiennych odpowiedzialnych za rozmiar"""
        self.max_size = max_size
        self.min_size = min_size
        # Ustawienia Entrez
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid):
        """Search for all records associated with a taxonomic ID."""
        print(f"Searching for records with taxID: {taxid}")
        try:
            # Najpierw pobierz informacje taksonomiczne
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            # Szukaj rekordów
            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)

            ##for record in search_results.keys():
                ##print(f"\t{record}") klu

                    #Count
                    #RetMax
                    #RetStart
                    #QueryKey
                    #WebEnv
                    #IdList
                    #TranslationSet
                    #TranslationStack
                    #QueryTranslation


            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name}")
                return None

            print(f"Found {count} records")

            # Zapisz wyniki wyszukiwania do późniejszego wykorzystania
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=10):
        """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            # Limit, aby zapobiec przeciążeniu serwera
            batch_size = min(max_records, 500)

            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )

            all_records = list(SeqIO.parse(handle, "genbank"))

            filtered_records = [
                record for record in all_records
                if self.min_size < len(record.seq) < self.max_size
            ]
            print(f"Filtered records: {len(filtered_records)}")
            # Surowy rekord GenBank
            # records_text = handle.read()

            return filtered_records

        except Exception as e:
            print(f"Error fetching records: {e}")
            return ""
    def generate_chart(self,records):
        pyplot.bar(records['accession'], records['sequence'],color='red',edgecolor='black')
        pyplot.ylabel('sequence length')
        pyplot.xlabel('accession')
        pyplot.savefig('wykres_sekwencji.pdf', format='pdf', bbox_inches='tight')
        pyplot.show()

def main():
    # Uzyskaj dane uwierzytelniające
    email = input("Enter your email address for NCBI: ")
    api_key = '113dfe4392781f84489d30d39a15522f7208' ##input("Enter your NCBI API key: ")
    #Ustaw długość sekwencji
    min = int(input("Enter minimum length of sequence to fetch: "))
    max = int(input("Enter maximum length of sequence to fetch: "))

    # Utwórz obiekt retriever
    retriever = NCBIRetriever(email, api_key,max,min)

    # Uzyskaj taxid od użytkownika
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")

    # Szukaj rekordów
    count = retriever.search_taxid(taxid)

    if not count:
        print("No records found. Exiting.")
        return

    # Pobierz kilka pierwszych rekordów jako próbkę
    print("\nFetching sample records...")
    sample_records = retriever.fetch_records(start=0, max_records=5)

    if sample_records:
        array = []
        for record in sample_records:
            array.append({
                "accession": record.id,
                "sequence": len(record.seq),
                "description": record.description
            })

    dateFrame = pd.DataFrame(array)
    dateFrame.to_csv('NCBI.csv', index=False)


    print(f"Saved sample records to NCBI.csv")

    ## Zapis danych do pdf
    retriever.generate_chart(records=dateFrame)


if __name__ == "__main__":
    main()
