# Config the directory

Move/Mount the directory which contains databases and ontologies files into accurate directory.

# Neo4j Settings

CAUTION: You need to comment `dbms.directories.import`, if you want to import files from anywhere in the filesystem.

If you run neo4j in a docker, you need to mount the related directory into docker.

```
# This setting constrains all `LOAD CSV` import files to be under the `import` directory. Remove or comment it out to
# allow files to be loaded from anywhere in the filesystem; this introduces possible security problems. See the
# `LOAD CSV` section of the manual for details.
#dbms.directories.import=/var/lib/neo4j/import
```
