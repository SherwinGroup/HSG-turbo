
my_list = ['a1', 'b1', 'c2', 'd2', 'e3', 'f3', 'g2', 'h1', 'i4', 'j4', 'k2']

new_list = []
for index in xrange(len(my_list)):
    try:
        temp = my_list.pop(0)
        print 'I just pooped', temp
    except:
        print 'I borked'
        break
    print 'the remaining list is ', my_list
    for spec in list(my_list):
        print 'bout to test ', spec
        if temp[-1] == spec[-1]:
            temp += spec
            my_list.remove(spec)
    new_list.append(temp)
    print 'after that round, the list is', my_list
print new_list