class myStr(object):
    def __init__(self, st):
        self.string = st
    def __getitem__(self, key):
        return self.string[key]
    def __str__(self):
        return self.string
    def __add__(self, other):
        return myStr(self.string + other.string)




my_list = ['a1', 'b1', 'A1', 'B2', 'c2', 'd2', 'e3', 'f3', 'g2', 'h1', 'i4', 'j4', 'k2']
my_list = ['a1', 'b1', 'A1', 'B2', 'c2', 'd2', 'e3', 'f3', 'i4', 'j4']
my_list = [myStr(i) for i in my_list]

#new_list = []
#for index in xrange(len(my_list)):
#    try:
#        temp = my_list.pop(0)
#        print 'I just pooped', temp
#    except:
#        print 'I borked'
#        break
#    print 'the remaining list is ', my_list
#    for spec in list(my_list):
#        print 'bout to test ', spec
#        if temp[-1] == spec[-1]:
#            temp += spec
#            my_list.remove(spec)
#    new_list.append(temp)
#    print 'after that round, the list is', my_list
#print new_list

new_list = []
for index in xrange(len(my_list)):
    try:
        temp = my_list.pop(0)
    except:
        break
    print "temp is: {}".format(temp)
    for spec in my_list:
        print "\tspec is: {}".format(spec)
        if temp[-1] == spec[-1]:
            temp += spec
            print "\t\tadded"
            my_list.remove(spec)
    new_list.append(temp)
#    print 'after that round, the list is', my_list
print new_list