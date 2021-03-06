import os
import logging

#- subset of desiutil.log.get_logger, to avoid desiutil dependency

_loggers = dict()
def get_logger(level=None, path=None, timestamps=False):
    '''Returns a logger, which can be used like logger.info('some message') etc.

    INPUTS: level ... None --> use environment variable
                      otherwise string like 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'FATAL', 'CRITICAL'

            path  ... None --> do not log to disk
                      otherwise string of where to save all log messages (identical console printouts will also happen)

            timestamps ... boolean, to turn on inclusion of timestamps in log messages

    OUTPUTS: logger instance
    '''
    if level is None:
        level = os.getenv('DESI_LOGLEVEL', 'INFO').upper()

    if level == 'DEBUG':
        loglevel = logging.DEBUG
    elif level == 'INFO':
        loglevel = logging.INFO
    elif level == 'WARN' or level == 'WARNING':
        loglevel = logging.WARNING
    elif level == 'ERROR':
        loglevel = logging.ERROR
    elif level == 'FATAL' or level == 'CRITICAL':
        loglevel = logging.CRITICAL
    else:
        raise ValueError('Unknown log level {}; should be DEBUG/INFO/WARNING/ERROR/CRITICAL'.format(level))

    if level not in _loggers:
        logger = logging.getLogger('desimeter.'+level)
        logger.setLevel(loglevel)

        #- handler and formatter code adapted from
        #- https://docs.python.org/3/howto/logging.html#configuring-logging

        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(loglevel)

        # optionally create file handler, similarly
        if path:
            fh = logging.FileHandler(filename=path, mode='a', encoding='utf-8')
            fh.setLevel(loglevel)
        else:
            fh = None

        # create formatter
        kwargs = {'fmt': '%(levelname)s:%(filename)s:%(lineno)s:%(funcName)s:%(message)s'}
        if timestamps:
            kwargs['fmt'] = '%(asctime)s:' + kwargs['fmt']
            kwargs['datefmt'] = '%Y%m%dT%H%M%S%z'
        formatter = logging.Formatter(**kwargs)

        # add formatter to ch
        ch.setFormatter(formatter)
        if fh:
            fh.setFormatter(formatter)

        # add handlers to logger
        logger.addHandler(ch)
        if fh:
            logger.addHandler(fh)

        _loggers[level] = logger

    return _loggers[level]
