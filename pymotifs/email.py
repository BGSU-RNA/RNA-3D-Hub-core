from mailer import Mailer
from mailer import Message


def body(filename):
    if not filename:
        return 'No log file produced'

    bad_lines = []
    with open(filename, 'rb') as raw:
        for line in raw:
            if line.startswith('WARNING') or line.startswith('ERROR') or \
                    line.startswith('CRITICAL'):
                bad_lines.append(line)

    return '\n'.join(bad_lines)


def send(stage, log_file=None, error=None):
    """Send an email based upon the log file.

    If no log file is provided or it is False then the email's body will not
    include the log lines to examine, and it will not be attached. If an error
    object is given this will be used to indicate the stage failed.

    :stage: Name of the stage that was run.
    :log_file: Filename for the log file
    :error: Exception object that occurred.
    """

    config = stage.config['email']

    status = 'Succeeded'
    if error:
        status = 'Failed'

    msg = Message(From=config['from'], To=config['to'])
    msg.Subject = config['subject'].format(stage=stage, status=status)
    msg.Body = body(log_file)
    if log_file:
        msg.attach(log_file)

    mail = Mailer(config['server']['name'])
    mail.send(msg)
