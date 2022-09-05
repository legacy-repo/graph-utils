import pandas as pd
import re
import logging
import os.path
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser

logger = logging.getLogger(__name__)


class UniProtParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'UniProt'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def build_stats(self):
        logger.info("Config file(%s): %s" % (self.config_fpath, self.config))
        stats = set()
        relationships_header = self.config['relationships_header']
        logger.info("Download release notes")
        self.parse_release_notes()
        logger.info("Download and parse fasta file...")
        stats = self.parse_fasta()

        # Proteins
        logger.info("Parse id mapping file...")
        stats.update(self.parse_idmapping_file())
        # Peptides
        logger.info("Parse uniprot peptides...")
        entities, relationships = self.parse_uniprot_peptides()
        entities_header = self.config['peptides_header']
        output_file = os.path.join(self.import_directory, "Peptide.tsv")
        logger.info("Generate file %s" % output_file)
        stats.update(self.print_single_file(entities, entities_header, output_file,
                     "entity", "Peptide", is_first=True, updated_on=self.updated_on))
        logger.info("Generate multiple files about relationships...")
        stats.update(self.print_multiple_relationships_files(
            relationships, relationships_header, self.import_directory, is_first=True, updated_on=self.updated_on))

        # Variants
        logger.info("Parse uniprot variants...")
        stats.update(self.parse_uniprot_variants())

        # Gene ontology annotation
        logger.info("Parse uniprot annotations for gene ontology...")
        relationships = self.parse_uniprot_annotations()
        relationships_header = self.config['go_header']
        logger.info("Generate multiple files about relationships...")
        stats.update(self.print_multiple_relationships_files(
            relationships, relationships_header, self.import_directory, is_first=True, updated_on=self.updated_on))

        # directory = os.path.join(self.database_directory, "UniProt")
        # self.remove_directory(directory)

        return stats

    def parse_release_notes(self):
        release_notes_url = self.config['release_notes']
        directory = os.path.join(self.database_directory, "UniProt")
        self.check_directory(directory)
        if self.download:
            self.download_db(release_notes_url, directory)

    def parse_fasta(self):
        stats = set()
        url = self.config['uniprot_fasta_file']
        entities_output_file = os.path.join(
            self.import_directory, "Amino_acid_sequence.tsv")
        rel_output_file = os.path.join(
            self.import_directory, "Protein_HAS_Sequence_Amino_acid_sequence.tsv")

        directory = os.path.join(self.database_directory, "UniProt")
        self.check_directory(directory)
        file_name = os.path.join(directory, url.split('/')[-1])

        if self.download:
            self.download_db(url, directory)

        ff = self.read_gzipped_file(file_name)
        records = self._parse_fasta(ff)
        num_entities = 0
        with open(entities_output_file, 'w', encoding='utf-8') as ef:
            logger.info("Generate entity file %s" % entities_output_file)
            ef.write('ID\theader\tsequence\tsize\tsource\n')
            with open(rel_output_file, 'w', encoding='utf-8') as rf:
                logger.info("Generate relationship file %s" % rel_output_file)
                rf.write('START_ID\tEND_ID\tTYPE\tsource\n')
                for i, batch in enumerate(self.batch_iterator(records, 1000)):
                    for record in batch:
                        identifier = record.id.split('|')[1]
                        header = record.id
                        sequence = str(record.seq)
                        sequence_len = len(str(sequence))
                        ef.write(identifier+"\t"+header+'\t'+sequence +
                                 '\t'+str(sequence_len)+'\tUniProt\n')
                        rf.write(identifier+'\t'+identifier +
                                 '\tHAS_SEQUENCE\tUniProt\n')
                        num_entities += 1

        stats.add(self._build_stats(num_entities, "entity",
                  "Amino_acid_sequence", "UniProt", entities_output_file, self.updated_on))
        stats.add(self._build_stats(num_entities, "relationships",
                  "HAS_SEQUENCE", "UniProt", rel_output_file, self.updated_on))

        return stats

    def parse_idmapping_file(self):
        regex_transcript = r"(-\d+$)"
        taxids = self.config['species']

        proteins_output_file = os.path.join(self.import_directory, "Protein.tsv")
        pdbs_output_file = os.path.join(
            self.import_directory, "Protein_structures.tsv")
        proteins = {}

        url = self.config['uniprot_id_url']
        directory = os.path.join(self.database_directory, "UniProt")
        self.check_directory(directory)
        file_name = os.path.join(directory, url.split('/')[-1])
        mapping_file = os.path.join(directory, 'mapping.tsv')
        if self.download:
            self.download_db(url, directory)

        fields = self.config['uniprot_ids']
        synonymFields = self.config['uniprot_synonyms']

        identifier = None
        transcripts = set()
        skip = set()
        is_first = True
        uf = self.read_gzipped_file(file_name)
        aux = {}
        stats = set()
        self.reset_mapping()
        with open(mapping_file, 'w', encoding='utf-8') as out:
            logger.info("Generate file %s" % proteins_output_file)
            logger.info("Generate file %s" % pdbs_output_file)
            logger.info("Generate multiple files about relationships...")
            for line in uf:
                data = line.rstrip("\r\n").split("\t")
                iid = data[0]
                if 'UniParc' in data:
                    if re.search(regex_transcript, iid):
                        transcripts.add(iid)
                    continue
                # data = line.decode('utf-8').rstrip("\r\n").split("\t")
                # iid = data[0]
                field = data[1]
                alias = data[2]

                if iid not in skip:
                    skip = set()
                    if re.search(regex_transcript, iid):
                        transcripts.add(iid)
                    if iid not in aux and iid.split('-')[0] not in aux:
                        if identifier is not None:
                            prot_info["synonyms"] = synonyms
                            aux[identifier].update(prot_info)
                            if "UniProtKB-ID" in aux[identifier] and "NCBI_TaxID" in aux[identifier]:
                                for synonym in synonyms:
                                    out.write(identifier+"\t"+synonym+"\n")
                                proteins[identifier] = aux[identifier]
                                for t in transcripts:
                                    proteins[t] = aux[identifier]
                                    for synonym in synonyms:
                                        out.write(t+"\t"+synonym+"\n")
                                if len(transcripts) > 0:
                                    proteins[identifier].update(
                                        {"isoforms": transcripts})
                                    transcripts = set()
                                aux.pop(identifier, None)
                                if len(proteins) >= 1000:
                                    entities, relationships, pdb_entities = self.format_output(
                                        proteins)
                                    stats.update(self.print_single_file(
                                        entities, self.config['proteins_header'], proteins_output_file, "entity", "Protein", is_first, self.updated_on))
                                    stats.update(self.print_single_file(
                                        pdb_entities, self.config['pdb_header'], pdbs_output_file, "entity", "Protein_structure", is_first, self.updated_on))
                                    stats.update(self.print_multiple_relationships_files(
                                        relationships, self.config['relationships_header'], self.import_directory, is_first, self.updated_on))
                                    is_first = False
                                    proteins = {}
                        identifier = iid
                        transcripts = set()
                        aux[identifier] = {}
                        prot_info = {}
                        synonyms = []
                    if field in fields:
                        if field == 'NCBI_TaxID':
                            if int(alias) not in taxids:
                                skip.add(identifier)
                                aux.pop(identifier, None)
                                identifier = None
                        if field in synonymFields:
                            synonyms.append(alias)
                        prot_info.setdefault(field, [])
                        prot_info[field].append(alias)

            uf.close()

        if len(proteins) > 0:
            entities, relationships, pdb_entities = self.format_output(proteins)
            stats.update(self.print_single_file(
                entities, self.config['proteins_header'], proteins_output_file, "entity", "Protein", is_first, self.updated_on))
            stats.update(self.print_single_file(
                pdb_entities, self.config['pdb_header'], pdbs_output_file, "entity", "Protein_structure", is_first, self.updated_on))
            stats.update(self.print_multiple_relationships_files(
                relationships, self.config['relationships_header'], self.import_directory, is_first, self.updated_on))

        self.mark_complete_mapping()

        return stats

    def format_output(self, proteins):
        entities = set()
        relationships = defaultdict(set)
        pdb_entities = set()
        for protein in proteins:
            accession = ""
            name = ""
            synonyms = []
            description = ""
            taxid = None
            pdb = []
            if "UniProtKB-ID" in proteins[protein]:
                accession = proteins[protein]["UniProtKB-ID"][0]
            if "Gene_Name" in proteins[protein]:
                name = proteins[protein]["Gene_Name"][0]
                for i in proteins[protein]["Gene_Name"]:
                    relationships[("Protein", "GENE_TRANSLATED_INTO")].add(
                        (i, protein, "GENE_TRANSLATED_INTO", 'UniProt'))
            if "RefSeq" in proteins[protein]:
                for i in proteins[protein]["RefSeq"]:
                    relationships[("Protein", "TRANSCRIPT_TRANSLATED_INTO")].add(
                        (i, protein, "TRANSCRIPT_TRANSLATED_INTO", 'UniProt'))
            if "synonyms" in proteins[protein]:
                synonyms = proteins[protein]["synonyms"]
            if "NCBI_TaxID" in proteins[protein]:
                taxid = proteins[protein]["NCBI_TaxID"][0]
                relationships[("Protein", "BELONGS_TO_TAXONOMY")].add(
                    (protein, int(taxid), "BELONGS_TO_TAXONOMY", 'UniProt'))
            if 'PDB' in proteins[protein]:
                pdb = proteins[protein]['PDB']
                for i in pdb:
                    pdb_entities.add(
                        (i, 'Protein_structure', 'Uniprot', 'http://www.rcsb.org/structure/{}'.format(i)))
                    relationships[("Protein", "HAS_STRUCTURE")].add(
                        (protein, i, "HAS_STRUCTURE", 'UniProt'))
            if "description" in proteins[protein]:
                description = proteins[protein]["description"]
            if "isoforms" in proteins[protein]:
                for i in proteins[protein]['isoforms']:
                    relationships[('Transcript', 'IS_ISOFORM')].add(
                        (i, protein, 'IS_ISOFORM', 'UniProt'))

            entities.add((protein, "Protein", accession, name,
                         ",".join(synonyms), description, int(taxid)))

        return entities, relationships, pdb_entities

    def print_single_file(self, data, header, output_file, data_type, data_object, is_first, updated_on):
        stats = set()
        df = pd.DataFrame(list(data), columns=header)
        stats.add(self._build_stats(len(data), data_type,
                  data_object, "UniProt", output_file, updated_on))
        with open(output_file, 'a', encoding='utf-8') as ef:
            df.to_csv(path_or_buf=ef, sep='\t',
                      header=is_first, index=False, quotechar='"',
                      line_terminator='\n', escapechar='\\')

        return stats

    def print_multiple_relationships_files(self, data, header, output_dir, is_first, updated_on):
        stats = set()
        for entity, relationship in data:
            df = pd.DataFrame(
                list(data[(entity, relationship)]), columns=header)
            output_file = os.path.join(
                output_dir, entity+"_"+relationship.lower() + ".tsv")
            stats.add(self._build_stats(len(data[(entity, relationship)]), 
                                        'relationships', relationship, "UniProt", 
                                        output_file, updated_on))
            with open(output_file, 'a', encoding='utf-8') as ef:
                df.to_csv(path_or_buf=ef, sep='\t',
                          header=is_first, index=False, quotechar='"',
                          line_terminator='\n', escapechar='\\')

        return stats

    def add_uniprot_texts(self, textsFile, proteins):
        with open(textsFile, 'r', encoding='utf-8') as tf:
            for line in tf:
                data = line.rstrip("\r\n").split("\t")
                protein = data[0]
                name = data[1]
                function = data[3]

                if protein in proteins:
                    proteins[protein].update({"description": function})

    def parse_uniprot_variants(self):
        data = defaultdict()
        variant_regex = r"(g\.\w+>\w)"
        chromosome_regex = r"(\w+)[p|q]"
        url = self.config['uniprot_variant_file']
        aa = self.config['amino_acids']
        entities = set()
        relationships = defaultdict(set)
        directory = os.path.join(self.database_directory, "UniProt")
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)
            
        known_variants = os.path.join(self.import_directory, "Known_variant.tsv")
        logger.info("Generate file %s" % known_variants)
        logger.info("Generate multiple files about relationships...")

        vf = self.read_gzipped_file(fileName)
        din = False
        i = 0
        stats = set()
        is_first = True
        for line in vf:
            line = line
            if not line.startswith('#') and not din:
                continue
            elif i <= 2:
                din = True
                i += 1
                continue
            data = line.rstrip("\r\n").split("\t")
            if len(data) > 9:
                gene = data[0]
                protein = data[1]
                pvariant = data[2]
                externalID = data[3]
                impact = data[4]
                clin_relevance = data[5]
                disease = data[6]
                chromosome_coord = data[8]
                original_source = data[13]
                ref = pvariant[2:5]
                pos = pvariant[5:-3]
                alt = pvariant[-3:]
                var_matches = re.search(variant_regex, data[9])
                chr_matches = re.search(chromosome_regex, chromosome_coord)
                if var_matches and chr_matches:
                    chromosome = 'chr'+chr_matches.group(1)
                    ident = chromosome+":"+var_matches.group(1)
                    altName = [externalID, data[5], pvariant, chromosome_coord]
                    if ref in aa and alt in aa:
                        altName.append(aa[ref]+pos+aa[alt])
                    pvariant = protein+"_"+pvariant
                    entities.add((ident, "Known_variant", pvariant, externalID, ",".join(
                        altName), impact, clin_relevance, disease, original_source, "UniProt"))
                    if chromosome != 'chr-':
                        relationships[('Chromosome', 'known_variant_found_in_chromosome')].add(
                            (ident, chromosome.replace('chr', ''), "VARIANT_FOUND_IN_CHROMOSOME", "UniProt"))
                    if gene != "":
                        relationships[('Gene', 'known_variant_found_in_gene')].add(
                            (ident, gene, "VARIANT_FOUND_IN_GENE", "UniProt"))
                    if protein != "":
                        relationships[('Protein', 'known_variant_found_in_protein')].add(
                            (ident, protein, "VARIANT_FOUND_IN_PROTEIN", "UniProt"))

                    if len(entities) >= 1000:
                        stats.update(self.print_single_file(entities, self.config['variants_header'],
                                                            known_variants, "entity", "Known_variant", is_first, self.updated_on))
                        stats.update(self.print_multiple_relationships_files(relationships, self.config['relationships_header'], self.import_directory, is_first, self.updated_on))
                        entities = set()
                        relationships = defaultdict(set)

                        is_first = False

        if len(entities) > 0:
            stats.update(self.print_single_file(entities, self.config['variants_header'], known_variants,
                                                "entity", "Known_variant", is_first, self.updated_on))
            stats.update(self.print_multiple_relationships_files(
                relationships, self.config['relationships_header'], self.import_directory, is_first, self.updated_on))
            del (entities)
            del (relationships)

        return stats

    def parse_uniprot_annotations(self):
        roots = {'F': 'Molecular_function',
                 'C': 'Cellular_component', 'P': 'Biological_process'}
        url = self.config['uniprot_go_annotations']
        relationships = defaultdict(set)
        directory = os.path.join(self.database_directory, "UniProt")
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)

        af = self.read_gzipped_file(fileName)
        for line in af:
            line = line
            if line.startswith('!'):
                continue
            data = line.rstrip("\r\n").split("\t")
            identifier = data[1]
            go = data[4]
            evidence = data[6]
            root = data[8]
            if root in roots:
                root = roots[root]
                relationships[(root, 'associated_with')].add(
                    (identifier, go, "ASSOCIATED_WITH", evidence, 5, "UniProt"))

        return relationships

    def parse_uniprot_peptides(self):
        file_urls = self.config['uniprot_peptides_files']
        entities = set()
        relationships = defaultdict(set)
        directory = os.path.join(self.database_directory, "UniProt")
        self.check_directory(directory)
        for url in file_urls:
            fileName = os.path.join(directory, url.split('/')[-1])
            if self.download:
                self.download_db(url, directory)
            first = True
            with open(fileName, 'r', encoding='utf-8') as f:
                for line in f:
                    if first:
                        first = False
                        continue

                    data = line.rstrip("\r\n").split("\t")
                    peptide = data[0]
                    accs = data[6].split(",")
                    is_unique = True
                    if len(accs) > 1:
                        is_unique = False
                    entities.add(
                        (peptide, "Peptide", "tryptic peptide", is_unique))
                    for protein in accs:
                        relationships[("Peptide", 'belongs_to_protein')].add(
                            (peptide, protein, "BELONGS_TO_PROTEIN", "UniProt"))
        return entities, relationships
