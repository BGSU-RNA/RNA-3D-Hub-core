from sender import Mail
from sender import Message
from sender import Attachement


def send(stage, error=None):
    config = stage.config['email']

    status = 'Succeeded'
    if error:
        status = 'Failed'

    bad_lines = []
    with open(log_file, 'rb') as raw:
        all_lines = []
        for line in raw:
            all_lines.append(line)
            if line.startswith('WARNING') or line.startswith('ERROR') or \
                    line.startswith('CRITICAL'):
                bad_lines.append(line)

        attachement = Attachement("output.log", "text/plain",
                                  '\n'.join(all_lines))

    msg = Message(config['subject'].format(stage=stage.name, status=status))
    msg.fromaddr = config['from']
    msg.body = '\n'.join(bad_lines)
    msg.to = config['to']
    msg.attach(attachement)

    mail = Mail(config['server']['name'], **config['server'].get('options'))
    mail.send(msg)
