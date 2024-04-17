from distutils.core import setup, Extension
import sysconfig

def main():
	CFLAGS = ['-g', '-Wall', '-std=c99', '-fopenmp', '-mavx', '-mfma', '-pthread', '-O3']
	LDFLAGS = ['-fopenmp']

	module = Extension(name = 'numc',
					extra_compile_args = CFLAGS,
					extra_link_args = LDFLAGS,
                    sources = ['numc.c', 'matrix.c'])

	# Use the setup function we imported and set up the modules.
	# You may find this reference helpful: https://docs.python.org/3.6/extending/building.html

	setup (name = 'numc',
       version = '1.0',
       description = 'This is CS 61c project 4 Numc Library',
       author = 'Annamira OToole and Teddi Worledge',
    #    author_email = 'annamira@berkeley.edu, tworledge@berkeley.edu',
       ext_modules = [module])

if __name__ == "__main__":
    main()

# library_dirs = ['/usr/local/lib'], include_dirs = ['/usr/local/include'], libraries = ['tcl83'],