project('stream2mmCIF', ['c'],
        version: '0.0.1',
        license: 'GPL3+',
        default_options: ['buildtype=debugoptimized'])


# Dependencies
cc = meson.get_compiler('c')
crystfel_dep = dependency('crystfel')
cif_dep = dependency('libcif')

# GUI
executable('stream2mmCIF',
           ['stream2mmCIF.c'],
           dependencies: [crystfel_dep, cif_dep],
           install: true)
