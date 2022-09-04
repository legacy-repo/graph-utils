import click
from builder.databases.databases_controller import database
from builder.ontologies.ontologies_controller import ontology

knowledge_graph = click.CommandCollection(sources=[database, ontology])