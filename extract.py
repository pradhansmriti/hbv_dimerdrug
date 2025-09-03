file='trajlammps.dat'
with open(file) as f:
    a =f.readlines()
indices=[i for i, j in enumerate(a) if 'LAMMPS' in j]
indices.append(len(a))
print(indices)
for i in range(3000,len(indices)-1):
    initial=indices[i]
    final=indices[i+1]
    with open('frame%s.dat'%i,'w') as fi:
        for index in range(initial,final):
            fi.write(a[index])
print('%d files'%(len(indices)-1))

