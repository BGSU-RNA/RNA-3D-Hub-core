"""This module contains some basic utilities for sending emails to the correct
users. It will build an email containing logging information and use that as
text for the email. This can probably be cleaned up and minimized with clever
use of better logging handlers. I think.
"""
import sys
import gzip
import shutil

from mailer import Mailer
from mailer import Message

from pymotifs import core

BODY = """
Input
-----
{input}

Argv
----
{args}

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
        """Compute the body text to use. This will contain infromation about
        the issues found during the run.

        Parameters
        ----------
        filename : str
            The path to the log file.

        Returns
        -------
        lines : str
            The configuration lines.
        """
        if not filename:
            return 'No log file produced'

        lines = []
        with open(filename, 'rb') as raw:
            for line in raw:
                if any(line.startswith(level) for level in self.levels):
                    lines.append(line)
        return ''.join(lines)

    def config_lines(self, config):
        """Create lines about the configuration used. The goal is to make it
        possible for some reading the email to recreate the input used.

        Parameters
        ----------
        config : dict
            The dictonary containing all options and configuration values.

        Returns
        -------
        lines : str
            The configuration lines.
        """
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
        id_lines = ' '.join('%s' % id for id in (ids or []))
        return BODY.format(options=config_lines,
                           issues=issues,
                           input=id_lines,
                           args=' '.join(sys.argv))

    def compressed_log(self, filename):
        """Create a compressed version of the given filename. This is useful
        for sending emails as it provides a compressed version of the file that
        can be attached to emails.

        Parameters
        ----------
        filename : str
            The filename to compress.

        Returns
        -------
        filename : str
            The compressed filename
        """
        if not filename:
            return None
        outname = filename + '.gz'
        with open(filename, 'rb') as raw, gzip.open(outname, 'wb') as out:
            shutil.copyfileobj(raw, out)
        return outname

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
            Body=self.body(name, log_file, **kwargs)
        )

        compressed_log = self.compressed_log(log_file)
        if compressed_log:
            msg.attach(compressed_log)

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
