print "Uploading the php files"

s='rsync -av '+conf.t.dir+'/ '+conf.g.username+',amtoolbox@web.sourceforge.net:/home/project-web/amtoolbox/htdocs/doc/'
print '   '+s
os.system(s)    
