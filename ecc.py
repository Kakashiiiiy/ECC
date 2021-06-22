class Curve(object):
	def __init__(self,p,a,b,x,y,q):
		super(Curve, self).__init__()
		self.p = p
		self.a = a
		self.b = b
		self.x = x
		self.y = y
		self.q = q

	def oncurve(self,x1,y1):
		x2 = (x1**3+ self.a *x1+ self.b)% self.p
		y2 = y1**2% self.p
		if x2 == y2:
			return True
		return False

class Point(object):
	def __init__(self, x,y, curve):
		super(Point, self).__init__()
		self.x = x
		self.y = y
		self.curve = curve

	def oncurve(self):
		return self.curve.oncurve(self.x,self.y)

	def pointdouble(self):
		if self.y == 0:
			return Point("O","O",self.curve)
		s = (3 * self.x**2 +self.curve.a)%self.curve.p
		s = s * inv_euklid(self.y*2,self.curve.p)
		x_new = (s**2 - 2*self.x) % self.curve.p
		y_new = (s*(self.x-x_new)-self.y) % self.curve.p
		return Point(x_new,y_new,self.curve)

	def pointaddition(self,p2):
		if self.y == -p2.y%self.curve.p and self.x == p2.x:
			return Point("O","O",self.curve)
		s = (p2.y-self.y) * inv_euklid(p2.x-self.x,self.curve.p)%self.curve.p
		x_new = (s**2 - self.x-p2.x) % self.curve.p
		y_new = (s*(self.x-x_new)-self.y) % self.curve.p
		return Point(x_new,y_new,self.curve)

	def return_coordinates(self):
		return self.x,self.y

	def return_inverse(self):
		return Point(self.x,-self.y,self.curve)

	def __add__(self,obj):
		if self.x == 'O':
			return obj
		if obj.x == 'O':
			return self
		if self==obj:
			ret = self.pointdouble()
		else:
			ret = self.pointaddition(obj)
		return ret
	def __str__(self):
		return "x: %s\ny: %s" % (hex(self.x), hex(self.y)) if type(self.x) is int else "x: O\ny: O"
	def __eq__(self,obj):
		return (self.x==obj.x%obj.curve.p and self.y==obj.y%obj.curve.p)
	def __sub__(self,obj):
		if obj.x == 'O':
			return self
		inv = Point(obj.x,-obj.y%obj.curve.p,obj.curve)
		return self.__add__(inv)
	def __mul__(self, a):
		return naf(self,a)
	def __rmul__(self, a):
		return self.__mul__(self,a)

def double_and_add(p,a):
	exp = bin(a)
	value = Point(p.x,p.y,p.curve)
	
	for i in range(3, len(exp)):
	    value = value + value 
	    if(exp[i:i+1]=='1'):
	        value = value + p
	return value

def calc_naf_representation(exponent):
	exp = bin(exponent)[2:]
	m = len(exp)
	loop_var = exponent
	naf_array = [0]*(m+1)
	naf_array[m]=2
	i=0
	for i in range(0,m+1):
		if loop_var == 0:
			break	
		x_ = loop_var%4
		if x_ == 1:
			naf_array[i] = 2 -x_ 
			loop_var = loop_var - 1
		elif x_ == 3:
			naf_array[i] = 2 - x_
			loop_var = loop_var + 1
		else:
			naf_array[i] = 0
		loop_var //= 2

	if naf_array[m]==2:
		naf_array = naf_array[:-1]
		m-= 1

	return naf_array


def naf(p,a):
	naf_array = calc_naf_representation(a)
	m = len(naf_array)-1
	value = Point(p.x,p.y,p.curve)
	m -=1

	while m >= 0:
		value = value+value
		if naf_array[m] == 1:
			value = value + p
		elif naf_array[m] == -1:
			value = value - p
		m-=1
	return value


def inv_euklid(a, b):
    x0, x1, y0, y1 = 0, 1, 1, 0
    a = a%b
    m = b
    while a != 0:
        q = b // a
        t = a
        a = b%a 
        b = t

        t = y0 
        y0 = y1 
        y1 = t - q * y1

        t = x0
        x0 = x1
        x1 = t - q * x1
    return x0%m

if __name__ == '__main__':
	
	curve = Curve(0xA9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377,0x7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9,0x26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6,0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262,0x547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997,0xA9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7)
	point = Point(0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262,0x547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997,curve)

	point_2 = Point(0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262,0x557c5fa5de13e4bea66dc47689226fa8abc4b110a73891d3c3f5f355f069e9e0,curve)
	#for i in range(10000):	
	inf = Point('O','O',curve)
	print(inverse*0xFAFAFAFAFAFAFAFAFAFAFAFAFAF%0xA9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7)
