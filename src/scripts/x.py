import string

def render( para_dic ): # , output_root):
  template_str = open('./template.html').read()

  for k, v in para_dic.items():
    print( str( k ) )
    template_str = template_str.replace( str('$'+k),  str(v) )
  print( template_str)

render( {'who': 'biodancer'})
