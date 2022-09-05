import os
import json
import logging
import logging.config as log_config
import click
import builder
from builder.databases import config as config_mod
from builder.databases.parsers import parsers
from joblib import Parallel, delayed

config_dir = os.path.dirname(os.path.abspath(builder.__file__))

# print(sys.path)
with open(os.path.join(config_dir, "log.config")) as f:
    config = json.load(f)
    log_config.dictConfig(config)

logger = logging.getLogger(__name__)


@click.group()
def database():
    pass


def parse_database(import_directory, database_directory, database, 
                   config_file=None, download=True, skip=True):
    stats = set()
    Parser = parsers.get(database, None)
    if Parser:
        parser = Parser(import_directory, database_directory, config_file, download, skip)
        stats = parser.build_stats()
    return stats


class NotSupportedAction(Exception):
    pass


@database.command(help="Make the databases and related graph files.")
@click.option('--db-dir', '-d', required=True,
              type=click.Path(exists=True, dir_okay=True),
              help="The directory which saved the downloaded database files.")
@click.option('--output-dir', '-o', required=True,
              type=click.Path(exists=True, dir_okay=True),
              help="The directory which saved the graph files.")
@click.option('--config', '-c', required=False,
              type=click.Path(exists=True, file_okay=True),
              help="The config file related with database.")
@click.option('--database', required=True, type=click.Choice(parsers.keys()),
              help="Which databases (you can specify the --database argument multiple times)?", multiple=True)
@click.option('--n-jobs', '-n', required=False,
              help="Hom many jobs?", default=4)
@click.option('--download/--no-download', default=False, help="Whether download the source file(s)?")
@click.option('--skip/--no-skip', default=True, help="Whether skip the existing file(s)?")
def make_database(output_dir, db_dir, database, config, download, n_jobs, skip):
    if config and len(database) > 1:
        raise NotSupportedAction("Cannot support a single config file with several databases.")

    all_databases = database
    valid_databases = list(
        filter(lambda database: database in parsers.keys(), all_databases))
    invalid_databases = list(
        filter(lambda database: database not in parsers.keys(), all_databases))
    if len(invalid_databases) > 0:
        logger.warn("%s databases (%s) is not valid, skip them.",
                    len(invalid_databases), invalid_databases)
    stats = Parallel(n_jobs=n_jobs)(delayed(parse_database)(
        output_dir, db_dir, database, config, download, skip) for database in valid_databases)
    allstats = {val if type(sublist) == set else sublist 
                for sublist in stats for val in sublist}
    logger.info("Stats: %s" % allstats)
    return allstats

@database.command(help="Print the default config file.")
@click.option('--database', required=True, type=click.Choice(parsers.keys()),
              help="Which database?", multiple=False)
def print_config(database):
    config_dir = os.path.dirname(os.path.abspath(config_mod.__file__))
    config_file = os.path.join(config_dir, "%s.yml" % database)
    with open(config_file, 'r') as f:
        print(f.read())
    

if __name__ == "__main__":
    main = click.CommandCollection(sources=[database])
    main()
