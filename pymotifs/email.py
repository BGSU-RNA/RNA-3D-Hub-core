from mailer import Mailer
from mailer import Message

from pymotifs import core

BODY = """
Input
-----
{input}

Options
-------
{options}

Issues
------
{issues}
"""


class Emailer(core.Base):
    """A class that handles the logic of sending an email. This will parse the
    log file to get the required lines, and send the mail if needed.
    """

    levels = set(['WARNING', 'ERROR', 'CRITICAL'])
    """Set of logging levels to put in email body."""

    skip_keys = set(['exclude', 'email', 'config', 'engine'])
    """Set of config keys to not put in email body."""

    def body_lines(self, filename):
        if not filename:
            return 'No log file produced'

        lines = []
        with open(filename, 'rb') as raw:
            for line in raw:
                if any(line.startswith(level) for level in self.levels):
                    lines.append(line)
        return ''.join(lines)

    def config_lines(self, config):
        lines = []
        for key, value in config.items():
            if key not in self.skip_keys:
                lines.append('{key}: {value}'.format(key=key, value=value))
        return '\n'.join(lines)

    def body(self, name, filename, ids=None, **kwargs):
        """Generate the body of the email.

        :param string filename: Name of the log file to create a body for.
        :returns: A string to use a body.
        """

        config_lines = self.config_lines(kwargs)
        issues = self.body_lines(filename)
        id_lines = ' '.join('%s' % id for id in ids)
        return BODY.format(options=config_lines, issues=issues, input=id_lines)

    def message(self, name, log_file=None, error=None, **kwargs):
        """Create an email for the given stage based upon the log file.

        If no log file is provided or it is False then the email's body will
        not include the log lines to examine, and it will not be attached. If
        an error object is given this will be used to indicate the stage
        failed.

        :stage: Name of the stage that was run.
        :log_file: Filename for the log file
        :error: Exception object that occurred.
        :returns: A Mailer message to send.
        """

        status = 'Succeeded'
        if error:
            status = 'Failed'

        kwargs = {'stage': name, 'status': status}
        subject = self.config['email']['subject'].format(**kwargs)
        to_address = kwargs.get('send_to', None) or self.config['email']['to']
        msg = Message(
            From=self.config['email']['from'],
            To=to_address,
            Subject=subject,
            Body=self.body(log_file, **kwargs)
        )
        if log_file:
            msg.attach(log_file)

        return msg

    def mailer(self, **kwargs):
        """Create the mailer used to send the email.

        :returns: The Mailer to use.
        """
        return Mailer(self.config['email']['server'])

    def __call__(self, name, **kwargs):
        """Send the email if needed. If this is a dry run no mail is sent and
        the email is just returned instead. If we have been configured to send
        no email then no email is built and None is returned. Otherwise the
        result of sending the email is returned.

        :name: The name of the stage that was run.
        """

        if self.config['email'].get('send', True) is False:
            return None

        mailer = self.mailer()
        msg = self.message(name, **kwargs)

        if kwargs.get('dry_run'):
            self.logger.info("Dry run so no mail sent")
            return msg

        return mailer.send(msg)
