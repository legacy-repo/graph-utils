import zipfile
import re
import logging
from os import path as os_path
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser

logger = logging.getLogger(__name__)


class CancerGenomeInterpreterParser(BaseParser):
    def __init__(self, import_directory, database_directory, download=True, skip=True) -> None:
        self.database_name = 'CancerGenomeInterpreter'
        config_dir = os_path.dirname(os_path.abspath(config.__file__))
        self.config_fpath = os_path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, download, skip)

    def parse(self):
        variant_regex = r"(\D\d+\D)$"
        regex = r"(chr\d+)\:g\.(\d+)(\w)>(\w)"
        url = self.config['cancerBiomarkers_url']
        entities_header = self.config['entities_header']
        relationships_headers = self.config['relationships_headers']
        amino_acids = self.config['amino_acids']
        mapping = self.get_mapping_from_ontology(ontology="Disease")
        drugmapping = self.get_mapping_for_entity("Drug")
        protein_mapping = self.get_multiple_mapping_for_entity("Protein")

        fileName = self.config['cancerBiomarkers_variant_file']
        relationships = defaultdict(set)
        entities = set()
        directory = os_path.join(
            self.database_directory, "CancerGenomeInterpreter")
        self.check_directory(directory)
        zipFile = os_path.join(directory, url.split('/')[-1])

        if self.download:
            self.download_db(url, directory)
        with zipfile.ZipFile(zipFile) as z:
            if fileName in z.namelist():
                with z.open(fileName, 'r') as responses:
                    first = True
                    for line in responses:
                        if first:
                            first = False
                            continue
                        data = line.decode('utf-8').rstrip("\r\n").split("\t")
                        gene_variant = data[0].split(':')
                        if len(gene_variant) < 2:
                            continue
                        gene = gene_variant[0]
                        variants = gene_variant[1].split(',')
                        #alterationType = data[1]
                        response = data[3]
                        drugs = data[10].split(';')
                        #status = data[11].split(';')
                        evidence = data[12]
                        tumors = data[16].split(';')
                        publications = data[17].split(';')
                        identifier = data[21]
                        prot_variant = data[22]
                        matches = re.match(regex, identifier)
                        alternative_names = [identifier]
                        if matches is not None:
                            cpra = matches.groups()
                            chromosome, position, reference, alternative = cpra
                            variant = chromosome+":g."+position+reference+">"+alternative
                            if prot_variant != "":
                                prot_variant = prot_variant.split(':')[1]
                                alternative_names.append(prot_variant)

                        valid_variants = []
                        if gene in protein_mapping:
                            for protein in protein_mapping[gene]:
                                for variant in variants:
                                    match = re.search(variant_regex, variant)
                                    if match:
                                        if variant[0] in amino_acids and variant[-1] in amino_acids:
                                            valid_variant = protein + '_p.' + \
                                                amino_acids[variant[0]] + ''.join(
                                                    variant[1:-1]) + amino_acids[variant[-1]]
                                            valid_variants.append(
                                                valid_variant)
                                            entities.add((valid_variant, "Clinically_relevant_variant",  ",".join(
                                                alternative_names), chromosome, position, reference, alternative, "", "", "CGI"))
                                            relationships["known_variant_is_clinically_relevant"].add(
                                                (valid_variant, valid_variant, "KNOWN_VARIANT_IS_CLINICALLY_RELEVANT", "CGI"))

                        for drug in drugs:
                            if drug.lower() in drugmapping:
                                drug = drugmapping[drug.lower()]
                            elif drug.split(" ")[0].lower() in drugmapping:
                                drug = drugmapping[drug.split(" ")[0].lower()]
                            elif " ".join(drug.split(" ")[1:]).lower() in drugmapping:
                                drug = drugmapping[" ".join(
                                    drug.split(" ")[1:]).lower()]
                            relationships["targets"].add(
                                (drug, gene, "CURATED_TARGETS", evidence, response, ",".join(tumors), "curated", "CGI"))

                            for valid_variant in valid_variants:
                                relationships["targets_clinically_relevant_variant"].add(
                                    (drug, valid_variant, "TARGETS_CLINICALLY_RELEVANT_VARIANT", evidence, response, "".join(tumors), "curated", "CGI"))

                        for tumor in tumors:
                            if tumor.lower() in mapping:
                                tumor = mapping[tumor.lower()]
                                for valid_variant in valid_variants:
                                    relationships["associated_with"].add(
                                        (valid_variant, tumor, "ASSOCIATED_WITH", "curated", "curated", "CGI", len(publications)))

        # self.remove_directory(directory)

        return (entities, relationships, entities_header, relationships_headers)

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header, relationships_headers = self.parse()
        entity_outputfile = os_path.join(
            self.import_directory, "cgi_clinically_relevant_variant.tsv")
        self.write_entities(
            entities, entities_header, entity_outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(self.database_name,
                    "clinically_relevant_variant", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "clinically_relevant_variant",
                                    self.database_name, entity_outputfile, self.updated_on))
        for relationship in relationships:
            cgi_outputfile = os_path.join(
                self.import_directory, "cgi_" + relationship + ".tsv")
            header = ['START_ID', 'END_ID', 'TYPE']
            if relationship in relationships_headers:
                header = relationships_headers[relationship]
            self.write_relationships(
                relationships[relationship], header, cgi_outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(
                self.database_name, relationship, len(relationships[relationship])))
            stats.add(self._build_stats(len(
                relationships[relationship]), "relationships", relationship, self.database_name, cgi_outputfile, self.updated_on))
        logger.info("Done Parsing database {}".format(self.database_name))
        return stats
