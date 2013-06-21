from fabric.api import cd
from fabric.api import env
from fabric.api import run
from fabric.api import task


env.hosts = ['anton@rna.bgsu.edu']
env.dir = '/Code/RNA-3D-Hub_dev'
env.branch = 'dev'


@task
def prod():
    env.dir = '/Code/RNA-3D-Hub'
    env.branch = 'master'


@task
def deploy():
    with cd(env.dir):
        run("git reset --hard")
        run("git pull origin %s" % env.branch)
