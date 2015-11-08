using UnityEngine;
using System.Collections.Generic;

//----------------------------------------------------------------------------------------------------------------------
///  SPHSolver class, this will be used to calculate the forces of a particle with the navier stokes equations


class SPHSolver
{

    /// @brief default ctor
    public SPHSolver()
    {
        m_smoothingLength = 1;
        m_externalForcePush = true;
    }

    //----------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------
    ///  Calculates the forces for a particle using leap frog method of integration
    ///  _currentParticle - The particle we are going to calculate the forces for
    ///  _particleIndex - Vector of all our particles if we want to brut force or the neighbour particles for SPH
    ///  _timeStep - the time step between updates
    public void calcForces(particle _currentParticle, List<particle> _particleIndex, float _timeStep)
    {
        float densityCurrent = 0;
        var gravitiy = new Vector3 (0f,-9.8f,0f);
        var pressureGrad = Vector3.zero;
        var viscosity = Vector3.zero;
        var viscosityVector = Vector3.zero;
        var Acceleration = Vector3.zero;
        for(int i = 0; i<_particleIndex.Count;i++){
        // Calculate desity of all particles
            if(_currentParticle!=_particleIndex[i]){
                densityCurrent += _particleIndex[i].getMass()*calcDensityKern(_currentParticle, _particleIndex[i]);
            }
        }
        _currentParticle.setDensity(densityCurrent);
        var pressTemp  = Vector3.zero;
        var visTemp  = Vector3.zero;
        float pressureCurrent, pressureNeighbour, pressure;

        for(int i=0; i<_particleIndex.Count;i++)
        {
            // Calcualate pressure
            pressureCurrent =   m_gasConstant  * (_currentParticle.getDensity() - _currentParticle.getRestDensity());
            pressureNeighbour = m_gasConstant  * (_particleIndex[i].getDensity() - _particleIndex[i].getRestDensity());
            pressure = ((pressureCurrent /(pressureCurrent*pressureCurrent)) + (pressureNeighbour /(pressureNeighbour*pressureNeighbour)))*(_particleIndex[i].getMass());
            pressTemp = calcPressureKern(_currentParticle,_particleIndex[i]);
            //pressTemp*= pressure; //!omfg not an assignment
            pressureGrad +=(pressTemp);
            // Calculate viscosiy vector
            viscosityVector = _particleIndex[i].getVel() -  _currentParticle.getVel();
            viscosityVector *=(_particleIndex[i].getMass());
            viscosityVector /=(_particleIndex[i].getDensity());
            // Calculate viscosiy
            visTemp = calcViscosityKern(_currentParticle,_particleIndex[i]);
            //visTemp = visTemp.cross(viscosityVector);   //!omfg not an assignment
            viscosity +=(visTemp);
        }
        viscosity *=m_viscosityCoefficient;
        //Calcualate external force if the is one   
        var Forces = Vector3.zero;
        if(m_externalForceStrenth != 0 && m_externalForceRadius != 0)
        {
            Forces = _currentParticle.getPos() - m_externalForcePos;
            if (Forces.length()<=m_externalForceRadius && Forces.length() > 0)
            {
                //find the direction the force is in
                Forces.Normalize();
                //Forces *=((m_externalForceRadius-Forces.length())/m_externalForceRadius);
                Forces *=(m_externalForceStrenth);
                if(m_externalForcePush == false)
                {

                    Forces =-Forces;
                }
            }
            else
            {
                Forces.set(0f,0f,0f);
            }
        }

        // Calculate our acceleration
        Acceleration = gravitiy - pressureGrad + viscosity + Forces;    //tweak center

        //---------------leap frog integration------------------------
        //Calculate velocity
        var VelHalfBack =   _currentParticle.getVel() -  _currentParticle.getAcc() * _timeStep/2;
        var VelHalfForward = VelHalfBack + Acceleration * _timeStep;
        var Velocity = (VelHalfBack + VelHalfForward) * 0.5f;
        //Calculate our new position
        var Position = _currentParticle.getPos() + VelHalfForward *(_timeStep);

        _currentParticle.setVel(Velocity);
        _currentParticle.setPos(Position);

//---------------Debuging----------------
//    Debug.Log("the viscosity is "+"["<<viscosity.m_x+","<<viscosity.m_y+","<<viscosity.m_z+"]");
//    Debug.Log("the pressure grad is "+"["<<pressureGrad.m_x<<","<<pressureGrad.m_y<<","<<pressureGrad.m_z<<"]");
//    Debug.Log("the accelleration is "+"["<<Acceleration.m_x<<","<<Acceleration.m_y<<","<<Acceleration.m_z<<"]");
//    Debug.Log("the velocity is "+"["<<Velocity.m_x+","<<Velocity.m_y<<","<<Velocity.m_z<<"]");

    }
    //----------------------------------------------------------------------------------------------------------------------
    ///  A mutatator for our smoothing length, h variable in navier stokes equations
    /// _smoothingLength - the smoothing length
    public void setSmoothingLength(float _smoothingLength)
    {
        m_smoothingLength = _smoothingLength;
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief returns our Smoothing Length if we need to query it
    public float getSmoothingLength()
    {
        return m_smoothingLength;
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Set our viscosity coefficient for our navier stokes equations
    public void setVisCoef(float _x)
    {
        m_viscosityCoefficient = _x;
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Initialises the density for our particles
    /// @param _particleIndex - The vector of all our particles we want to initialize the density for
    public void initDensity(List<particle> _particleIndex)
    {
        float densityCurrent = 0.0f;
        for(int i=0; i<_particleIndex.Count;i++)
        {
        // Calculate desity of all particles
            for(int j=0; j<_particleIndex.Count;j++){
                if(i!=j){
                    densityCurrent += _particleIndex[j].getMass()*calcDensityKern(_particleIndex[i], _particleIndex[j]);
                }
            }
            //We cant have 0 densities or it will effect out calculations this is fine for bruce mesh but for anything arbetory if could make problems
            //this isnt precise but for approximations it will do
            if(densityCurrent==0.0) densityCurrent = 1.0f;

            //Set the densties for our particle
            _particleIndex[i].setDensity(densityCurrent);
            _particleIndex[i].setRestDensity(densityCurrent);

            //reinit our temp density for next loop
            densityCurrent = 0.0f;
        }
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief set the gas constant for our water simulation
    /// @param _x - the value we wish to set our constant to
    public void setGasConstant(float _x)
    {
        m_gasConstant = _x;
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Add an external force to our SPH calculations. acts like wind
    /// @param _pos - The source position of the force
    /// @param _forceRadius - the area that the force has influence over. This has a linear falloff
    /// @param _forceStrength - How strong your want your force to be. This is just a scaler.
    /// @param _push - a bool so we know weather we are pushing the particles away or sucking them towards us
    public void setExternalForce(Vector3 _pos, float _forceRadius, float _forceStrength, bool _push)
    {
        m_externalForcePos = _pos;
        m_externalForceRadius = _forceRadius;
        m_externalForceStrenth = _forceStrength;
        m_externalForcePush = _push;
    }
    //----------------------------------------------------------------------------------------------------------------------

    /// @brief calculates the smoothing kernal for density, used in our SPH equations
    /// @param _currentParticle the particle we wish to test for
    /// @param _neighbour the particle we wish to test our _current particle against
    /// @return the weight of influence the particle has on our _currentParticle
    protected float calcDensityKern(particle _currentParticle, particle _neighbour)
    {
        Vector3 r = _currentParticle.getPos() - _neighbour.getPos();
        if(r.length() > m_smoothingLength)
        {
            return 0;
        }
        float densityKern = (315 / (64 * Mathf.PI * Mathf.Pow(m_smoothingLength, 9))) * Mathf.Pow(((m_smoothingLength * m_smoothingLength) - (r.length() * r.length())), 3);
        return densityKern;
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief calculates the smoothing kernal for density, used in our SPH equations
    /// @param _currentParticle the particle we wish to test for
    /// @param _neighbour the particle we wish to test our _current particle against
    /// @return the weighted pressued gradient the particle has on our _currentParticle
    protected Vector3 calcPressureKern(particle _currentParticle,particle _neighbour)
    {
        Vector3 r = _currentParticle.getPos() - _neighbour.getPos();
        if (r.length()>m_smoothingLength)
        {
            return Vector3.zero;
        }
        float pressureKern = -(945/(32*Mathf.PI * Mathf.Pow(m_smoothingLength,9))) * Mathf.Pow(((m_smoothingLength*m_smoothingLength) - (r.length()*r.length())),3);
        return r *(pressureKern);
    }

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief calculates the smoothing kernal for density, used in our SPH equations
    /// @param _currentParticle the particle we wish to test for
    /// @param _neighbour the particle we wish to test our _current particle against
    /// @return the weighted viscosity vector the particle has on our _currentParticle
    protected Vector3 calcViscosityKern(particle _currentParticle, particle _neighbour)
    {
        Vector3 r = _currentParticle.getPos() - _neighbour.getPos();

        if(r.length() > m_smoothingLength)
        {
            return Vector3.zero;
        }
        float viscosityKern = -(945/(32*Mathf.PI * Mathf.Pow(m_smoothingLength,9))) * Mathf.Pow(((m_smoothingLength*m_smoothingLength) - (r.length()*r.length())),3) * ((3*(m_smoothingLength*m_smoothingLength)) - 7*(m_smoothingLength*m_smoothingLength));
        return r * viscosityKern;
    }

    //----------------------------------------------------------------------------------------------------------------------
    
    private float m_smoothingLength;
    //----------------------------------------------------------------------------------------------------------------------
    
    private float m_viscosityCoefficient;
    //----------------------------------------------------------------------------------------------------------------------
    
    private float m_gasConstant;
    //----------------------------------------------------------------------------------------------------------------------
    
    private Vector3 m_externalForcePos;
    //----------------------------------------------------------------------------------------------------------------------
    
    private float m_externalForceRadius;
    //----------------------------------------------------------------------------------------------------------------------
   
    private float m_externalForceStrenth;
    //----------------------------------------------------------------------------------------------------------------------
    
    private bool m_externalForcePush;
    //----------------------------------------------------------------------------------------------------------------------

}

