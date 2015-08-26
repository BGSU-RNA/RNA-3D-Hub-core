from mailer import Mailer
from mailer import Message

from pymotifs import core


class Emailer(core.Base):
    """A class that handles the logic of sending an email. This will parse the
    log file to get the required lines, and send the mail if needed.
    """

    def body(self, filename):
        """Generate the body of the email.
        """

        if not filename:
            return 'No log file produced'

        bad_lines = []
        with open(filename, 'rb') as raw:
            for line in raw:
                if line.startswith('WARNING') or line.startswith('ERROR') or \
                        line.startswith('CRITICAL'):
                    bad_lines.append(line)

        return '\n'.join(bad_lines)

    def message(self, name, log_file=None, error=None, **kwargs):
        """Create an email for the given stage based upon the log file.

        If no log file is provided or it is False then the email's body will
        not include the log lines to examine, and it will not be attached. If
        an error object is given this will be used to indicate the stage
        failed.

        :stage: Name of the stage that was run.
        :log_file: Filename for the log file
        :error: Exception object that occurred.
        """

        status = 'Succeeded'
        if error:
            status = 'Failed'

        msg = Message(From=self.config['from'], To=self.config['to'])
        msg.Subject = self.config['subject'].format(stage=name, status=status)
        msg.Body = self.body(log_file)
        if log_file:
            msg.attach(log_file)

        return msg

    def mailer(self, **kwargs):
        """Create the mailer used to send the email.
        """

        return Mailer(self.config['server']['name'])

    def __call__(self, name, **kwargs):
        """Send the email if needed. If this is a dry run no mail is sent and
        the email is just returned instead. If we have been configured to send
        no email then no email is built and None is returned. Otherwise the
        result of sending the email is returned.

        :name: The name of the stage that was run.
        """

        if not self.config['email'].get('send'):
            self.logger.info("Configured to not send mail")
            return None

        msg = self.message(name, **kwargs)

        if kwargs.get('dry_run'):
            self.logger.info("Dry run so no mail sent")
            return msg

        mailer = self.mailer()
        return mailer.send(msg)
