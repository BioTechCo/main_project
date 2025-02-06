import sys
import logging


logging_config = dict(
    version=1,
    formatters={
        'verbose': {
            'format': ("[%(asctime)s] %(levelname)s "
                       "[%(name)s:%(lineno)s] %(message)s"),
            'datefmt': "%d/%b/%Y %H:%M:%S",
        },
        'simple': {
            'format': '%(levelname)s %(message)s',
        },
    },
    handlers={
        'file_handler': {
            'class': 'logging.handlers.RotatingFileHandler',
            'formatter': 'verbose',
            'level': logging.DEBUG,
            'filename': 'logs/app.log',
            'maxBytes': 52428800,
            'backupCount': 7,
        },
        'console': {
            'class': 'logging.StreamHandler',
            'level': 'DEBUG',
            'formatter': 'simple',
            'stream': sys.stdout,
        },
    },
    loggers={
        'train_helper': {
            'handlers': ['file_handler', 'console'],
            'level': logging.DEBUG,
        },
        'simple_model': {
            'handlers': ['file_handler', 'console'],
            'level': logging.DEBUG,
        },
        'process_norm': {
            'handlers': ['file_handler', 'console'],
            'level': logging.DEBUG,
        },
    }
)
