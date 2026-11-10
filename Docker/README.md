# Docker

The image installs conda packages from `environment-linux-64.lock`, a lockfile generated from `environment.yml`.

To change conda packages, edit `environment.yml`, then regenerate the lockfile from this directory:

```bash
conda-lock lock -f environment.yml -p linux-64 --kind explicit --filename-template 'environment-{platform}.lock'
```

Commit both files. The lockfile only covers linux-64.
