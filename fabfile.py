from fabric.api import cd
from fabric.api import env
from fabric.api import run
from fabric.api import task


env.hosts = ['pipeline@rnatest']
env.dir = '/usr/local/code/pipeline'
env.branch = 'dev'


@task
def prod():
    env.dir = ['pipeline@rna']
    env.branch = 'master'


@task
def deploy():
    with cd(env.dir):
        run("git reset --hard")
        run("git pull origin %s" % env.branch)
