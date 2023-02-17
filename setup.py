import sys, os
from glob import glob
from setuptools import setup,find_packages

snap_files = []
snapdir = './nasap/scripts/snap/**'

snap_directories = [f for f in glob(snapdir, recursive=True) if os.path.isdir( f )]
# print( snap_directories )
for directory in snap_directories:
  snap_files.append((directory, [os.path.join(directory, f) for f in os.listdir(directory) if not os.path.isdir(os.path.join(directory, f))]) )

gencore_files = ['nasap/scripts/gencore/gencore']
bg2bw_files = ['nasap/scripts/ucsc/bedGraphToBigWig']

data_files = [('nasap/scripts/', glob('nasap/scripts/[!_]*.py', recursive=True) ),
  ('nasap/scripts/', glob('nasap/scripts/[!_]*.bash', recursive=True) ),
  ('nasap/scripts/templates/', ['nasap/scripts/templates/template.html', 'nasap/scripts/templates/template_track.txt']),
  ('nasap/scripts/templates/static/', ['nasap/scripts/templates/static/vue.min.js']),
]
# print( snap_files )
data_files.extend(snap_files)
data_files.extend( gencore_files )
data_files.extend( bg2bw_files )

def get_version():
  with open( os.path.abspath(os.path.dirname(__file__) ) + '/nasap/__version__.py') as f:
    ver = f.read().split('version=')[1].replace('"', '').strip()
  return ver

def main():
  root = os.path.abspath(os.path.dirname(__file__))
  with open(os.path.join(root, 'requirements.txt')) as f:
    install_requires = f.read().splitlines()
  # linux depandance packages
  # linux_deps = ['fastp', 'bioawk', 'python', 'bowtie2', 'samtools', 'bedtools', 'deeptools']
  # for dep in linux_deps:
  #   is_dep = os.system(dep + ' --version')
  #   if is_dep != 0:
  #     apt install

  # if float(sys.version[:3])< 3.6:
  #   sys.stderr.write("CRITICAL: Python version must be >= 3.6x!\n")
  #   sys.exit(1)

  setup(
    name='nasap',
    version=get_version(),
    description='This is nASAP setup file',
    author='biodancer',
    author_email='szxszx@foxmail.com',
    url='https://github.com/biodancerwanghzi/nasap/',
    # packages= find_packages(exclude=["back"]),
    packages=['nasap'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    data_files= data_files,
    entry_points={
      'console_scripts': [
        'nasap=nasap.nasap:main',
        'batch_nasap=nasap.batch_nasap:main'
      ]
    },
    install_requires=install_requires,
    extras_require={
      'genomeBrowser': [
        'pyGenomeTrack'
    ]

    },
    # zip_safe = True
  )

if __name__ == '__main__':
  main()