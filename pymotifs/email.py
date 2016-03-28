import logging

from mailer import Mailer
from mailer import Message

from pymotifs import core

logger = logging.getLogger(__name__)


class Emailer(core.Base):
    """A class that handles the logic of sending an email. This will parse the
    log file to get the required lines, and send the mail if needed.
    """

    def body(self, filename):
        """Generate the body of the email.

        :param string filename: Name of the log file to create a body for.
        :returns: A string to use a body.
        """

        if not filename:
            return 'No log file produced'

        bad_lines = []
        with open(filename, 'rb') as raw:
            for line in raw:
                if line.startswith('WARNING') or line.startswith('ERROR') or \
                        line.startswith('CRITICAL'):
                    bad_lines.append(line)

        return ''.join(bad_lines)

    @property
    def mailer(self):
        """Create the mailer used to send the email.

        :returns: The Mailer to use.
        """
        return Mailer(self.config['email']['server'])

    def message(self, name, log_file=None, status=None):
        """Create an email for the given stage based upon the log file.

        If no log file is provided or it is False then the email's body will
        not include the log lines to examine, and it will not be attached. If
        an error object is given this will be used to indicate the stage
        failed.

        :param str stage: Name of the stage that was run.
        :param str log_file: Filename for the log file
        :param str status: The status to use in the subject line.
        :returns: A Mailer message to send.
        """

        kwargs = {'stage': name, 'status': status}
        subject = self.config['email']['subject'].format(**kwargs)
        msg = Message(
            From=self.config['email']['from'],
            To=self.config['email']['to'],
            Subject=subject,
            Body=self.body(log_file)
        )
        if log_file:
            msg.attach(log_file)

        return msg

    def success(self, name, log_file=None, dry_run=False):
        """Send a sucess message.
        """
        self.msg = self.message(name, status='Success', log_file=log_file)
        return self.__send__(self.msg)

    def fail(self, name, log_file=None, dry_run=False):
        """Send a failure message.
        """
        self.msg = self.message(name, status='Failed', log_file=log_file)
        return self.__send__(self.msg)

    def __send__(self, message, dry_run=False):
        """Send a a message. This takes some message and sends it. If an
        exception occurs None is returned and the exception is logged.
        """

        try:
            if dry_run:
                self.logger.info("Dry run so no mail sent")
            else:
                self.mailer.send(self.msg)
                return message
        except Exception as err:
            logger.exception(err)
            return None
