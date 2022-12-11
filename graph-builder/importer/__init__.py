import os
import re
import sys
import yaml
import click
import logging
import verboselogs
import coloredlogs

from urllib.parse import unquote
from importer import connector

verboselogs.install()
coloredlogs.install(
    fmt='%(asctime)s - %(module)s:%(lineno)d - %(levelname)s - %(message)s')
logger = logging.getLogger('importer')


def read_yaml(yaml_file):
    content = None
    with open(yaml_file, 'r') as stream:
        try:
            content = yaml.safe_load(stream)
        except yaml.YAMLError as err:
            raise yaml.YAMLError(
                "The yaml file {} could not be parsed. {}".format(yaml_file, err))
    return content


def connect_neo4j(db_url="localhost:7687", user="neo4j", password="password"):
    try:
        uri = "bolt://".format(db_url)
        driver = neo4j.GraphDatabase.driver(
            uri, auth=(user, password), encrypted=False)
    except Exception as err:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        sys_error = "{}, file: {},line: {}".format(
            sys.exc_info(), fname, exc_tb.tb_lineno)
        print("Database is not online")
        #raise Exception("Unexpected error:{}.\n{}".format(err, sys_error))

    return driver


def load_into_database(driver, queries, requester):
    """
    This function runs the queries provided in the graph database using a neo4j driver.

    :param driver: neo4j driver, which provides the connection to the neo4j graph database.
    :type driver: neo4j driver
    :param list[dict] queries: list of queries to be passed to the database.
    :param str requester: identifier of the query.
    """
    regex = r"file:\/\/\/(.+\.tsv)"
    result = None
    for query in queries:
        try:
            if "file" in query:
                matches = re.search(regex, query)
                if matches:
                    file_path = matches.group(1)
                    if os.path.isfile(unquote(file_path)):
                        result = connector.commitQuery(driver, query+";")
                        record = result.single()
                        if record is not None and 'c' in record:
                            counts = record['c']
                            if counts == 0:
                                logger.warning(
                                    "{} - No data was inserted in query: {}.\n results: {}".format(requester, query, counts))
                            else:
                                logger.info(
                                    "{} - Query: {}.\n results: {}".format(requester, query, counts))
                        else:
                            logger.info(
                                "{} - cypher query: {}".format(requester, query))
                    else:
                        logger.error(
                            "Error loading: File does not exist. Query: {}".format(query))
            else:
                result = connector.commitQuery(driver, query+";")
        except Exception as err:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            logger.error(
                "Loading: {}, file: {}, line: {} - query: {}".format(err, fname, exc_tb.tb_lineno, query))

    return result


@click.group()
def importer():
    pass


@importer.command(help="Import ontology files into graph database.")
@click.option('--ontology-dir', '-d', required=True,
              type=click.Path(exists=True, dir_okay=True),
              help="The directory which saved the parsed ontology files.")
@click.option('--db-url', '-d', required=True, help="Neo4j database url. Please contain database name, such as localhost:7687/default.")
@click.option('--db-username', '-u', default="neo4j", help="Neo4j database username.")
@click.option('--db-password', '-P', default="NeO4J", help="Neo4j database password.")
def ontology(ontology_dir, db_url, db_username, db_password):
    # TODO: how to notice user that the ontologies files need to be imported firstly.
    config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "config", "cypher.yml")
    importer_config = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "config", "importers.yml")
    queries = read_yaml(config_file)
    importers = read_yaml(importer_config)

    entities = [e.lower() for e in importers["ontology_entities"]]
    ontology_data_import_code = queries["IMPORT_ONTOLOGY_DATA"]['query']
    formated_codes = []
    for entity in entities:
        formated_codes.extend(ontology_data_import_code.replace(
            "ENTITY", entity.capitalize()).replace("IMPORTDIR", ontology_dir).split(";")[0:-1])
    mappings = importers["ontology_mappings"]
    mapping_import_code = queries["IMPORT_ONTOLOGY_MAPPING_DATA"]['query']
    for m in mappings:
        if m.lower() in entities:
            for r in mappings[m]:
                formated_codes.extend(mapping_import_code.replace("ENTITY1", m).replace(
                    "ENTITY2", r).replace("IMPORTDIR", ontology_dir).split(";")[0:-1])
    print("Dont loading ontologies.")

    driver = connect_neo4j(
        db_url=db_url, user=db_username, password=db_password)
    logger.info("Updating ontologies...")
    load_into_database(driver=driver, queries=formated_codes,
                       requester="ontologies")


ENTITY_MAP = {
    "Gene": ["IMPORT_GENE_DATA"],
    "Chromosome": ["IMPORT_CHROMOSOME_DATA"],
    "Transcript": ["IMPORT_TRANSCRIPT_DATA"],
    "Protein": ["IMPORT_PROTEIN_DATA", "IMPORT_PROTEIN_ANNOTATIONS", "IMPORT_MODIFIED_PROTEIN_ANNOTATIONS", "IMPORT_MODIFIED_PROTEINS"],
}

DATABASE_ENTITY_MAP = {
    "IMPORT_MODIFIED_PROTEIN_ANNOTATIONS": ["PhosphoSitePlus"],
    "IMPORT_MODIFIED_PROTEINS": ["SIGNOR", "PhosphoSitePlus"]
}


@importer.command(help="Import database files into graph database.")
@click.option('--database-dir', '-d', required=True,
              type=click.Path(exists=True, dir_okay=True),
              help="The directory which saved the parsed database files.")
@click.option('--which-entity', 'w', required=True, type=click.Choice(ENTITY_MAP.keys()), help="Which database you want to import.")
@click.option('--db-url', '-d', required=True, help="Neo4j database url. Please contain database name, such as localhost:7687/default.")
@click.option('--db-username', '-u', default="neo4j", help="Neo4j database username.")
@click.option('--db-password', '-P', default="NeO4J", help="Neo4j database password.")
def database(database_dir, which_entity, db_url, db_username, db_password):
    # TODO: how to notice user that the ontologies files need to be imported firstly.
    config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "config", "cypher.yml")
    queries = read_yaml(config_file)

    importers = ENTITY_MAP.get(which_entity)
    formated_codes = []
    for importer in importers:
        databases = DATABASE_ENTITY_MAP.get(importer)
        import_code = queries[importer]['query']
        # Maybe have several databases which provide the related files.
        if databases:
            for database in databases:
                formated_codes.extend(import_code.replace("IMPORTDIR", database_dir).replace(
                    "RESOURCE", database).split(";")[0:-1])
        else:
            formated_codes.extend(import_code.replace(
                "IMPORTDIR", database_dir).split(";")[0:-1])

    print("Dont loading %s's files." % which_entity)

    driver = connect_neo4j(
        db_url=db_url, user=db_username, password=db_password)
    logger.info("Updating %s..." % which_entity)
    load_into_database(driver=driver, queries=formated_codes,
                       requester=which_entity)


if __name__ == "__main__":
    importer()
