import subprocess

build = subprocess.run("python setup.py sdist bdist_wheel" , shell=True, capture_output=True, text=True)

print(build.stdout)
print(build.stderr)

publish = subprocess.run("twine upload dist/* --repository qchem" , shell=True, capture_output=True, text=True)

print(publish.stdout)
print(publish.stderr)