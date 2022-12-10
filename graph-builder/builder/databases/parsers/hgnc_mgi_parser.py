import os.path
import logging
import pandas as pd
import re
import verboselogs
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')

NA_MARKER = ""
# NA_MARKER = "NA"


class HGNC_MGI_Parser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'HGNC_MGI'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def parse(self):
        hgnc = self.parse_hgnc()
        mgi = self.parse_mgi()
        return [*hgnc[0], *mgi[0]], hgnc[1]

    def parse_hgnc(self):
        url = self.config['hgnc_url']
        entities = set()
        directory = os.path.join(self.database_directory, "HGNC_MGI")
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        taxid = 9606  # Homo sapiens
        entities_header = self.config['header']

        if self.download:
            self.download_db(url, directory)

        with open(fileName, 'r', encoding="utf-8") as df:
            first = True
            for line in df:
                if first:
                    first = False
                    continue
                data = line.rstrip("\r\n").split("\t")
                gene_symbol = data[1]
                gene_name = data[2]
                status = data[5]
                gene_family = data[12]
                synonyms = data[18:23]
                entrez_id = ",".join(list(filter(lambda item: re.match(r"[0-9]+", str(item)), synonyms))) or NA_MARKER
                ensembl_id = ",".join(list(filter(lambda item: re.match(r"EN.*", str(item)), synonyms))) or NA_MARKER
                # transcript = data[23]
                if status != "Approved":
                    continue

                entities.add((gene_symbol, "Gene", gene_name,
                             gene_family, entrez_id, ensembl_id, ",".join(synonyms), taxid))
                #relationships.add((geneSymbol, transcript, "TRANSCRIBED_INTO"))

        # self.remove_directory(directory)

        return entities, entities_header

    def merge_duplicated(self, filepath):
        df = pd.read_csv(filepath, delimiter="\t", dtype=object, index_col=False)
        renamed_d = df.rename(columns={
            "3. marker symbol": "gene_symbol",
            "4. marker name": "name",
            "5. genome build": "build",
            "6. Entrez gene id": "entrez_id",
            "7. NCBI gene chromosome": "chr",
            "8. NCBI gene start": "start",
            "9. NCBI gene end": "end",
            "10. NCBI gene strand": "strand",
            "11. Ensembl gene id": "ensembl_id",
        })

        selected_df = renamed_d[['gene_symbol', 'name', 'build',
                                 'entrez_id', 'chr', 'start', 'end', 'strand', 'ensembl_id']]
        return selected_df

        # def concat_func(x):
        #     return pd.Series({
        #         'name': ','.join(list(map(str, x['name']))),
        #         'build': ','.join(list(map(str, x['build']))),
        #         'entrez_id': ','.join(list(map(str, x['entrez_id']))),
        #         'chr': ','.join(list(map(str, x['chr']))),
        #         'start': ','.join(list(map(str, x['start']))),
        #         'end': ','.join(list(map(str, x['end']))),
        #         'strand': ','.join(list(map(str, x['strand']))),
        #         'ensembl_id': ','.join(list(map(str, x['ensembl_id']))),
        #     })

        # selected_df = renamed_d[['gene_symbol', 'name', 'build',
        #                          'entrez_id', 'chr', 'start', 'end', 'strand', 'ensembl_id']]
        # result = selected_df.groupby(selected_df['gene_symbol']).apply(
        #     concat_func).reset_index()
        # return result

    def parse_mgi(self):
        url = self.config['mgi_url']
        entities = set()
        directory = os.path.join(self.database_directory, "HGNC_MGI")
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        taxid = 10090  # Mus musculus
        entities_header = self.config['header']

        if self.download:
            self.download_db(url, directory)

        df = self.merge_duplicated(fileName)
        for (index, row) in df.iterrows():
            gene_symbol = row["gene_symbol"]
            gene_name = row["name"]
            gene_family = ""
            entrez_id = row["entrez_id"] if str(row["entrez_id"]) != "nan" else NA_MARKER
            ensembl_id = row["ensembl_id"] if str(row["ensembl_id"]) != "nan" else NA_MARKER
            other_synonyms = ""

            entities.add((gene_symbol, "Gene", gene_name,
                          gene_family, entrez_id, ensembl_id, other_synonyms, taxid))
            #relationships.add((geneSymbol, transcript, "TRANSCRIBED_INTO"))

        # self.remove_directory(directory)

        return entities, entities_header

    def build_stats(self):
        stats = set()
        entities, header = self.parse()
        outputfile = os.path.join(self.import_directory, "Gene.tsv")
        self.write_entities(entities, header, outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(
            self.database_name, "Gene", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "Gene",
                  self.database_name, outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
