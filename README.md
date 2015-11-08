# SPHsimntroduction

I created a simulation using Smooth particle hydrodynamics to create fluid simulation. I would apply external forces to the particles to see what the effects of the particles are.


Implementation

 The value of a force at a given position is interpolated from a discrete set of points. SPH derives from the integral interpolation below.


W(r.h) is the smoothing function
The SPH simulation is a collection of particles i that interact with another set of particles N, within radius of h of i, dictated by the smoothing function. Each particle i has a position r, a velocity, and a density p. 
   
protected Vector3 calcPressureKern(particle _currentParticle,particle _neighbour)
{  Vector3 r = _currentParticle.getPos() - _neighbour.getPos();
        if (r.length()>m_smoothingLength)
        {return Vector3.zero
float pressureKern = -(945/(32*Mathf.PI * Mathf.Pow(m_smoothingLength,9))) * Mathf.Pow(((m_smoothingLength*m_smoothingLength) - (r.length()*r.length())),3);
return r *(pressureKern);

protected Vector3 calcPressureKern(particle _currentParticle,particle _neighbour)
    {
        Vector3 r = _currentParticle.getPos() - _neighbour.getPos();
        if (r.length()>m_smoothingLength)
        { return Vector3.zero;}
        float pressureKern = -(945/(32*Mathf.PI * Mathf.Pow(m_smoothingLength,9))) * Mathf.Pow(((m_smoothingLength*m_smoothingLength) - (r.length()*r.length())),3);
        return r *(pressureKern)}
The steps to calculate SPH is to (1)calculate density at every particle position (2) calculate pressure forces, viscosity forces for each particle, and (3) update the acceleration structures. 
 
Spatial hashing
In SPH particles interact only with those particles that are within the circle of radius h. Most of the particles do not interact and interaction only occurs within the circle of radius. To optimise which are the most important collisions spatial hashing is used, positions in a scene are hashed to create a key, which specifies its position in a table. The purpose of this technique is so that positions close to each other are hashed to the same bucket in the hash table.  So in order to find neighbours of a particle at a certain position they can be looked up at its hashed position in the table and in the same bucket there should lie its neighbours. We would ideally like to take to take a position of the particle and convert it to a unique hash number. The hash number must be unique so as not to have two particles at different locations have the same positional cell.

protected int hashFunction(particle _particle)
    {
        int x_x = Mathf.FloorToInt(_particle.getPos().x / m_smoothingLength);
        int x_y = Mathf.FloorToInt(_particle.getPos().y / m_smoothingLength);
        int x_z = Mathf.FloorToInt(_particle.getPos().z / m_smoothingLength);
        return (x_x * m_primes[0] ^ x_y * m_primes[1] ^ x_z * m_primes[2]) & m_hashTableSize;
    }
This is the hash function that generates a key based on the particle position, the prime numbers are so that 





Results
showing the destruction of particle man

https://www.youtube.com/watch?v=BVkVxfogrFs&feature=youtu.be


