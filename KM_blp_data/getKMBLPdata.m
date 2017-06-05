datadir = '/Users/shoshievass/Dropbox/WorkSpaces/Discrete_Choice_Models/KM_blp_data';
cd(datadir) 


dlmwrite('cdid.txt',cdid,'delimiter','\t')
dlmwrite('price.txt',price,'delimiter','\t')
dlmwrite('x1.txt',x1,'delimiter','\t')
dlmwrite('x2.txt',x2,'delimiter','\t')
dlmwrite('IV1.txt',IV1,'delimiter','\t')
dlmwrite('id.txt',id,'delimiter','\t')
dlmwrite('firmid.txt',firmid,'delimiter','\t')
dlmwrite('share.txt',share,'delimiter','\t')
dlmwrite('outshr.txt',outshr,'delimiter','\t')
