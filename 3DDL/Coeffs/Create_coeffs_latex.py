""" This python script aims to create a latex file where all the tension condition coefficients are written """
import re

path_X = '/home/jordan/Documents/Memoire_figure_git/3DDL/coeff_X.txt'
path_Y = '/home/jordan/Documents/Memoire_figure_git/3DDL/coeff_Y.txt'
path_Z = '/home/jordan/Documents/Memoire_figure_git/3DDL/coeff_Z.txt'
path_K = '/home/jordan/Documents/Memoire_figure_git/3DDL/coeff_K.txt'



def coeff_array(Path):
    coeff_array = []
    with open(Path,'r') as coeff:
        for line in coeff.readlines():
            coeffs = line.split(',    ')
            for i in range(len(coeffs)):
                coeffs[i] = re.sub("""\\\,""", '', coeffs[i])
            coeffs.pop()
            coeff_array.append(coeffs)

    coeff_array = list(map(list, zip(*coeff_array)))

    return coeff_array

def get_var_name(*variable):
    '''gets string of variable name
    inputs
        variable (str)
    returns
        string
    '''
    if len(variable) != 1:
        raise Exception('len of variables inputed must be 1')
    try:
        return [k for k, v in locals().items() if v is variable[0]][0]
    except:
        return [k for k, v in globals().items() if v is variable[0]][0]

def create_list_tuple(l):
    l_name = get_var_name(l)
    ret_list = []
    for i in range(len(l)):
        ret_list.append(('{}_{{{}{}}}'.format(l_name[0],l_name[1],i+1),l[i]))
    return ret_list



coeff_X = coeff_array(path_X)
coeff_Y = coeff_array(path_Y)
coeff_Z = coeff_array(path_Z)
coeff_K = coeff_array(path_K)

lt ={'ay':0, 'az':0, 'ak':0, 'bx':0, 'bz':0, 'bk':0, 'cx':0, 'cy':0, 'ck':0, 'dx':0, 'dk':0, 'ey':0, 'ek':0, 'hz':0, 'hk':0}

bx = coeff_X[1]
lt['bx'] = create_list_tuple(bx)
cx = coeff_X[2]
lt['cx'] = create_list_tuple(cx)
dx = coeff_X[3]
lt['dx'] = create_list_tuple(dx)

ay = coeff_Y[0]
lt['ay'] = create_list_tuple(ay)
cy = coeff_Y[2]
lt['cy'] = create_list_tuple(cy)
ey = coeff_Y[4]
lt['ey'] = create_list_tuple(ey)

az = coeff_Z[0]
lt['az'] = create_list_tuple(az)
bz = coeff_Z[1]
lt['bz'] = create_list_tuple(bz)
hz = coeff_Z[5]
lt['hz'] = create_list_tuple(hz)

ak = coeff_K[0]
lt['ak'] = create_list_tuple(ak)
bk = coeff_K[1]
lt['bk'] = create_list_tuple(bk)
ck = coeff_K[2]
lt['ck'] = create_list_tuple(ck)
dk = coeff_K[3]
lt['dk'] = create_list_tuple(dk)
ek = coeff_K[4]
lt['ek'] = create_list_tuple(ek)
hk = coeff_K[5]
lt['hk'] = create_list_tuple(hk)

lts = {}
for item in lt:
    lts[item] = ''
    i=1
    for val in lt[item]:
        if i ==2 or (i==1 and item=='ck'):
            lts[item] += '{}={},\\\\ '.format(val[0], val[1])
            i=1
        else:
            lts[item] += '{}={}, '.format(val[0], val[1])
            i+=1
    lts[item] = lts[item][:-4]
#lts[item] = lts[item][:-4]
with open('/home/jordan/Documents/Memoire_jordan/tex/coeffs.tex','w') as tf:
      tf.write(
          '\\begin{align*} \n '+
          lts['ay']+'\\\\'+ lts['az'] + ',\\\\' + lts['ak'] + ',\\\\' + lts['bx'] + ',\\\\' + lts['bz'] +
          ',\\\\' + lts['bk']+ ',\\\\' + lts['cx']+ ',\\\\' + lts['cy']+ ',\\\\' + lts['ck']+ ',\\\\' + lts['dx']+
          ',\\\\' + lts['dk'] + ',\\\\' + lts['ey']+ ',\\\\' + lts['ek']+ ',\\\\' + lts['hz']+ ',\\\\' + lts['hk']+
          '\n\\end{align*}')

tf.close()
