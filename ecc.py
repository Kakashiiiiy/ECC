ka = 0x20A5B20E076E77984380CB49173F6ED7FDED87E645747133F63888907245E5D8 
kb = 0x63690612179A5742A7DB7003F0545E866CAF9DE086BF272A0E1827165381B399

#1.1 und 1.2 kann man sich schenken, da es dafür keine gesonderten Funktionen in Python braucht

class Curve(object):
	"""docstring for Curve"""
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
	"""docstring for Point"""
	def __init__(self, x,y, curve):
		super(Point, self).__init__()
		self.x = x
		self.y = y
		self.curve = curve

	def oncurve(self):
		return self.curve.oncurve(self.x,self.y)

	def pointdouble(self):
		s = (3 * self.x**2 +self.curve.a)%self.curve.p
		s = s * inv_euler(self.y*2,self.curve.p)
		x_new = (s**2 - 2*self.x) % self.curve.p
		y_new = (s*(self.x-x_new)-self.y) % self.curve.p
		return Point(x_new,y_new,self.curve)

	def pointadd(self,p2):
		s = (p2.y-self.y) * inv_euler(p2.x-self.x,self.curve.p)%self.curve.p
		x_new = (s**2 - self.x-p2.x) % self.curve.p
		y_new = (s*(self.x-x_new)-self.y) % self.curve.p
		return Point(x_new,y_new,self.curve)

	def return_coordinates(self):
		return self.x,self.y

	def return_inverse(self):
		return Point(self.x,-self.y,self.curve)

	def __add__(self,obj):
		if self==obj:
			ret = self.pointdouble()
		else:
			ret = self.pointadd(obj)
		return ret
	def __str__(self):
		return "x: %s\ny: %s" % (hex(self.x), hex(self.y))
	def __eq__(self,obj):
		return (self.x==obj.x and self.y==obj.y)
	def __sub__(self,obj):
		inv = Point(obj.x,-obj.y%obj.curve.p,obj.curve)
		return self.__add__(inv)
	def __mul__(self, a):
		return naf(self,a)

def double_add(p,a):
	exp = bin(a)
	value = Point(p.x,p.y,p.curve)
	
	for i in range(3, len(exp)):
	    value = value + value 
	    if(exp[i:i+1]=='1'):
	        value = value + p
	return value

def naf(p,a):
	exp = bin(a)
	value = Point(p.x,p.y,p.curve)
	m = len(exp)
	i = 2
	while exp[i%m]=='1':
		value = value + value
		i+=1
	value = value + p.return_inverse()
	#print(value)

	while i < m:
		value = value + value
		if exp[(i+1)%m]=='1':
			if exp[(i+2)%m]=='1':
				value = value+p
				i+=1
				while exp[i%m]=='1':
					value = value+value
					i+=1
				value = value + p.return_inverse()
				continue
		if exp[i%m]=='1':
			value = value + p
		i+=1
	return value

def exp_func(x, y, m): #x**y % m
    exp = bin(y)
    value = x
 
    for i in range(3, len(exp)):
        value = value * value % m
        if(exp[i:i+1]=='1'):
            value = value*x % m
    return value 

def inv_euklid(a, b):
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        q, b, a = b // a, a, b % a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return x0 

def inv_euler(a,b):
	return exp_func(a,b-2,b)

curve = Curve(0xA9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377,0x7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9,0x26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6,0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262,0x547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997,0xA9FB57DBA1EEA9BC3E660A909D838D718C397AA3B561A6F7901E0E82974856A7)
point = Point(0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262,0x547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997,curve)

print("doubleandadd:euler")
for i in range(10000):
	p123 = double_add(point,(ka*kb%curve.q))

