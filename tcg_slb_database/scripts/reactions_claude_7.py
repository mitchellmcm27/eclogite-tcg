
an = ['Feldspar', 'Anorthite']    
ab = ['Feldspar', 'Albite']
fs = ['Orthopyroxene', 'Ferrosilite']
oEn = ['Orthopyroxene', 'Enstatite']
mgts = ['Orthopyroxene', 'MgTschermaks']
oDi = ['Orthopyroxene', 'OrthoDiopside']
hed = ['Clinopyroxene', 'Hedenbergite']
di = ['Clinopyroxene', 'Diopside']
jd = ['Clinopyroxene', 'Jadeite']
cats = ['Clinopyroxene', 'CaTschermaks']
cEn = ['Clinopyroxene', 'Clinoenstatite']
alm = ['Garnet', 'Almandine']
pyp = ['Garnet', 'Pyrope']
gs = ['Garnet', 'Grossular']
qz = ['Quartz', 'Quartz']
ky = ['Kyanite', 'Kyanite']
mgsp = ['Spinel', 'MgSpinel']
herc = ['Spinel', 'Hercynite']
fo = ['Olivine', 'Forsterite']
fa = ['Olivine', 'Fayalite']

reactions_claude = [
    ([cEn], [oEn]),
    ([di], [oDi]),
    ([fs, oDi], [oEn, hed]),
    ([fs, mgts], [oEn, alm]),
    ([oDi, mgts], [oEn, cats]),
    ([oDi, mgts], [oEn, gs]),
    ([hed, cEn], [fs, di]),
    ([cEn, cats], [mgts, di]),
    ([cEn, cats], [di, pyp]),
    ([cEn, gs], [di, cats]),
    ([mgts, herc], [ky, fs, mgsp]),
    ([oDi, alm], [oEn, hed, gs]),
    ([mgts, hed], [oEn, cats, alm]),
    ([mgts, mgsp, fa], [oEn, herc]),
    ([di, alm], [oEn, hed, gs]),
    ([hed, pyp], [oEn, cats, alm]),
    ([oDi, alm], [fs, cEn, gs]),
    ([oDi, alm], [fs, cats, pyp]),
    ([mgts, hed], [fs, di, gs]),
    ([mgts, hed], [fs, cEn, gs]),
    ([mgts, hed], [fs, cats, pyp]),
    ([mgts, herc, fo], [fs, mgsp]),
    ([mgts, fa], [fs, herc, fo]),
    ([hed, pyp], [fs, di, gs]),
    ([di, alm], [fs, cEn, gs]),
    ([di, alm], [fs, cats, pyp]),
    ([hed, pyp], [fs, cEn, gs]),
    ([mgts, hed], [oDi, cats, alm]),
    ([hed, cEn, cats], [di, alm]),
    ([cats, mgsp, fa], [di, herc]),
    ([cats, herc, fo], [hed, mgsp]),
    ([cEn, cats, herc], [qz, hed, mgsp]),
    ([cEn, cats, herc], [ky, hed, mgsp]),
    ([cEn, cats, fa], [ky, hed, fo]),
    ([ab, oDi, mgts], [an, oEn, jd]),
    ([ab, mgts, di], [an, cEn, jd]),
    ([ab, di, pyp], [an, cEn, jd]),
    ([ab, cEn, gs], [an, di, jd]),
    ([ab, di, mgsp], [an, cEn, jd]),
    ([oDi, mgts, herc], [an, fs, mgsp]),
    ([cEn, cats, herc], [an, hed, mgsp]),
    ([oDi, mgts, herc], [fs, cats, mgsp]),
    ([oDi, mgts, fa], [fs, cats, fo]),
    ([oDi, mgts, herc], [fs, gs, mgsp]),
    ([oDi, mgts, fa], [fs, gs, fo]),
    ([cEn, cats, herc], [mgts, hed, mgsp]),
    ([cEn, cats, fa], [mgts, hed, fo]),
    ([cEn, cats, herc], [di, alm, mgsp]),
    ([cEn, cats, fa], [di, alm, fo]),
    ([cEn, cats, herc], [hed, pyp, mgsp]),
    ([cEn, cats, fa], [hed, pyp, fo]),
    ([cEn, gs, herc], [hed, cats, mgsp]),
    ([cEn, gs, fa], [hed, cats, fo]),
    ([cEn, cats, herc], [hed, alm, mgsp]),
    ([cEn, cats, fa], [hed, alm, fo]),
    ([ab, mgts, hed], [an, oEn, jd, alm]),
    ([ab, mgts, hed], [an, oEn, jd, fa]),
    ([ab, hed, pyp], [an, oEn, jd, alm]),
    ([ab, hed, mgsp], [an, oEn, jd, herc]),
    ([ab, oDi, alm], [an, fs, jd, pyp]),
    ([ab, oDi, herc], [an, fs, jd, mgsp]),
    ([ab, mgts, hed], [an, fs, jd, pyp]),
    ([ab, mgts, hed], [an, fs, jd, fo]),
    ([ab, di, alm], [an, fs, jd, pyp]),
    ([ab, di, herc], [an, fs, jd, mgsp]),
    ([ab, mgts, hed], [an, oDi, jd, alm]),
    ([ab, di, alm], [an, hed, cEn, jd]),
    ([ab, di, herc], [an, hed, cEn, jd]),
    ([ab, di, herc], [an, jd, alm, mgsp]),
    ([ab, di, herc], [an, jd, mgsp, fa]),
    ([ab, hed, mgsp], [an, jd, pyp, herc]),
    ([ab, hed, mgsp], [an, jd, herc, fo]),
    ([cEn, cats, jd, herc], [ab, hed, mgsp]),
]