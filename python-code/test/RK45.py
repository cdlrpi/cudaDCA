import numpy as np
import strtDCA as stDCA


def RK45_strtDCA(state,tspan,n,bs,js):

	h=tspan[1]
	q=np.zeros((n))
	qq=np.zeros((n))
	qqdot=np.zeros((n))
	qdot=np.zeros((n))
	q1=np.zeros((n))
	qdot1=np.zeros((n))
	q2=np.zeros((n))
	qdot2=np.zeros((n))
	q3=np.zeros((n))
	qdot3=np.zeros((n))
	q4=np.zeros((n))
	qdot4=np.zeros((n))
	qddot1=np.zeros((n))
	qddot2=np.zeros((n))
	qddot3=np.zeros((n))	
	qddot4=np.zeros((n))

	qq[:]=state[:n]
	qqdot[:]=state[n:]
	Y=np.zeros((2*n))
	#print(h)

	q1[:]=qq[:]
	qdot1[:]=qqdot[:]
	state[:n]=q1[:]
	state[n:]=qdot1[:]
	#print('state')
	#print(state+np.array((-np.pi/2,0,0,0)))
	Y1= stDCA.DCA(bs,js,state)
	qddot1=Y1[n:]
	#print("Y1")
	#print(Y1)

	q2[:]=qq+(1/2)*h*qdot1
	qdot2[:] = qqdot+ (1/2)*h*qddot1
	state[:n]=q2[:]
	state[n:]=qdot2[:]
	#print('state2')
	#print(state+np.array((-np.pi/2,0,0,0)))
	Y2=stDCA.DCA(bs,js,state)
	qddot2[:]=Y2[n:]
	#print("Y2")
	#print(Y2)

	#print("q")
	#print(qq)
	#print("qdot2")
	#print(qdot2)
	#print("qdot")
	#print(qqdot)
	#print("qddot2")
	#print(qddot2)
	q3[:]=qq+(1/2)*h*qdot2
	qdot3[:]=qqdot+(1/2)*h*qddot2
	
	state[:n]=q3[:]
	state[n:]=qdot3[:]
	#print('state3')
	#print(state+np.array((-np.pi/2,0,0,0)))
	Y3=stDCA.DCA(bs,js,state)
	qddot3[:]=Y3[n:]
	#print("Y3")
	#print(Y3)
	q4[:] = qq+h*qdot3
	qdot4[:]=qqdot+h*qddot3
	state[:n]=q4
	state[n:]=qdot4
	#print('state4')
	#print(state+np.array((-np.pi/2,0,0,0)))
	Y4=stDCA.DCA(bs,js,state)
	qddot4[:]=Y4[n:]

	#print("Y4")
	#print(Y4)
	q[:]=qq+(h/6)*(qdot1+2*qdot2+2*qdot3+qdot4)
	qdot[:]=qqdot+(h/6)*(qddot1+2*qddot2+2*qddot3+qddot4)
	Y[:n]=q
	Y[n:]=qdot
	state[:n]=qq[:]
	state[n:]=qqdot[:]
	#print("Y")
	#print(Y-np.array((np.pi/2,0,0,0)))
	return Y

