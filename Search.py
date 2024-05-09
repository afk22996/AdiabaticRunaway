import numpy as np

'''Modified binary search which returns a tuple of the two indices between which the target value should lie between. If the value is in the list, the returned tuple contains duplicate values
with both being the index of the target value. If the value is outside of the bounds of the list, a tuple containing either positive or negative infinities is used where the other value is unknown
to commmunicate that uncertainty'''
def binSearch(vals, left, right, target):
	while left < right:
		middle = left + (right - left)//2
		if(middle >= len(vals)-1):
			break

		if(vals[middle] == target): #target value is in array
			return (middle, middle)

		elif(vals[middle +1] > target and vals[middle] < target): #correct guess with the next value being the upper bound
			return (middle, middle+1)

		elif(vals[middle-1] < target and vals[middle] > target): #correct value with the previous value being the lower bound
			return (middle-1, middle)

		else:
			if(vals[middle] < target): #guess was too low
				left = middle + 1

			else: #guess was too high
				right = middle - 1

	if(right > len(vals)-1): #target was outside of upper bound of array
		return(len(vals)-1, np.infty)

	else: #target was outside of lower bound of array
		return(-np.infty,0)
		
if __name__ == '__main__':
	arr = [0,1,2,3,4,6,7]
	print(binSearch(arr, 0, len(arr), 7))