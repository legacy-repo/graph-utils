import os
import sys
import neo4j


def commitQuery(driver, query, parameters={}):
    result = None
    try:
        with driver.session() as session:
            result = session.run(query, parameters)
    except neo4j.exceptions.ClientError as err:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        sys_error = "{}, file: {},line: {}".format(
            sys.exc_info(), fname, exc_tb.tb_lineno)
        print("Connection error:{}.\n{}".format(err, sys_error))
    except Exception as err:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        sys_error = "{}, file: {},line: {}".format(
            sys.exc_info(), fname, exc_tb.tb_lineno)
        raise Exception("Connection error:{}.\n{}".format(err, sys_error))

    return result
