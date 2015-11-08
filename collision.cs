using UnityEngine;
using System.Collections.Generic;


//----------------------------------------------------------------------------------------------------------------------
/// @file collision
/// @brief This is my collision class, It will be used to test collisions of particles and update accordingly.


//----------------------------------------------------------------------------------------------------------------------

public class collision
{
    ///-------------------------------------------------------------------------------------------------
    /// @brief A functioin to run through all of the walls and test a selected particle for collision
    /// @param _testparticle The particle that we wish to test
    /// @param _timestep used in calculating the friction against the wall
    public void testCollisionWithWalls(particle _testParticle, float _timeStep)
    {
        //for all our walls calculate collision
        for(int i = 0; i < m_walls.Count; i++)
        {
            calculateWallCollision(_testParticle,_timeStep,m_walls[i].plane);
        }
    }

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Add a wall to our scene
    /// @param _normal the direction we wish for our wall to point in
    /// @param _center the postion we wish for our wall to be
    /// @param _draw a bool for if we want to draw
    /// @param _drawSizeX how big we want to draw our wall in X scale
    /// @param _drawSizeY how big we want to draw our wall in Y scale
    /// @param _drawSizeZ how big we want to draw our wall in Z scale
    public void addWall(Vector3 _normal, Vector3 _center, bool _draw, float _drawSizeX, float _drawSizeY, float _drawSizeZ)
    {
        //render part
        if (m_prefab == null)
        {
            m_prefab = Resources.Load("models/wall") as GameObject;
        }

        GameObject m_go = GameObject.Instantiate(m_prefab, _center, Quaternion.LookRotation(_normal)) as GameObject;

        m_go.name = "Wall" + m_gos.Count;
        m_go.transform.localScale = new Vector3(_drawSizeX, _drawSizeY, _drawSizeZ);
        m_gos.Add(m_go); 
        ////////////////////////////////////////////////////////////
        /// collision logical part
        /// 
        wall tempWall = new wall();
        tempWall.draw = _draw;
        if (_draw)
        {
            tempWall.drawSizeX = _drawSizeX;
            tempWall.drawSizeY = _drawSizeY;
            tempWall.drawSizeZ = _drawSizeZ;
        }

        tempWall.center = _center;
        tempWall.plane = new Vector4(_normal.x,_normal.y,_normal.z,_center.x*_normal.x + _center.y*_normal.y + _center.z*_normal.z);
        m_walls.Add(tempWall);
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief If particles are inside each other moves them apart
    /// @param _testParticle the particle we want to test
    /// @param _particles our array of particles we wish to test against
    public void initParticleCollision(particle _testParticle, List<particle> _particles)
    {
        Vector3  newPosP, newPosQ, currentPosQ,currentPosP,N;

        float radius = _testParticle.getRadius();
        for (int i=0; i<_particles.Count;i++)
        {
            N = _testParticle.getPos() - _particles[i].getPos();
            if (N.sqrMagnitude > 0 && N.sqrMagnitude < radius)
            {
                currentPosP = _testParticle.getPos();
                currentPosQ = _particles[i].getPos();
                newPosP = (currentPosP  + (N  * (0.5f)));
                newPosQ = (currentPosQ  - (N  * (0.5f)));
                _testParticle.setPos(newPosP);
                _particles[i].setPos(newPosQ);
            }
        }
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Set the coefficient of restitution for our walls
    /// @param _x Value we wish to set to
    public void setCoefOfRest(float _x)
    {
        m_coefOfRest = _x;
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Set the coefficient of friction for our walls
    /// @param _x Value we wish to set to
    public void setCoefOfFric(float _x)
    {
        m_coefOfFric = _x;
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief draws our walls to the scene
    /// @param _trans Our global transformation stack of our scene so we can retrieve model matrix
    /// @param _cam The camera of our scene to retrieve view and perspective matrix
    public void drawWalls(Transform tr, Camera _cam)
    {
        for(int i = 0; i < m_gos.Count; i++)
        {
            if(m_walls[i].draw && !m_gos[i].activeSelf)
            {
                m_gos[i].SetActive(true);
            }
            else if(!m_walls[i].draw && m_gos[i].activeSelf)
            {
                m_gos[i].SetActive(false);
            }
        }
    }

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A mutator to change the position of our walls
    /// @param _wall the wall in our array we wish to change
    /// @param _normal the normal we wish to change it to
    /// @param _center the center of the wall we wish to change to
    public void setWall(int _wall, Vector3 _normal, Vector3 _center, bool _draw, float _drawSizeX , float _drawSizeY , float _drawSizeZ)
    {
        if (_wall>= m_walls.Count)
        {
            Debug.LogError("wall not in range");
            return;
        }
        if(_normal.sqrMagnitude > 0)
        {
            m_walls[_wall].Center(_center);
            m_walls[_wall].Draw(_draw);
            m_walls[_wall].Plane(new Vector4(_normal.x,_normal.y,_normal.z,_center.x*_normal.x + _center.y*_normal.y + _center.z*_normal.z));
            m_walls[_wall].DrawSizeX(_drawSizeX);
            m_walls[_wall].DrawSizeY(_drawSizeY);
            m_walls[_wall].DrawSizeZ(_drawSizeZ);
        }
        else
        {
            m_walls[_wall].Center(_center);
            m_walls[_wall].Draw(_draw);
            m_walls[_wall].Plane(new Vector4(m_walls[_wall].plane.x, m_walls[_wall].plane.y, m_walls[_wall].plane.z, _center.x * _normal.x + _center.y * _normal.y + _center.z * _normal.z));
            m_walls[_wall].DrawSizeX(_drawSizeX);
            m_walls[_wall].DrawSizeY(_drawSizeY);
            m_walls[_wall].DrawSizeZ(_drawSizeZ);
        }
    }


    /// @brief calculates and updates a particles position and velocity depending on if it has a collision with the ground
    /// @param _testParticle the particle we wish to update
    /// @param _timeStep the time step between our calcuations, this is used in calculting the friction
    protected void calculateWallCollision(particle _testParticle, float _timeStep, Vector4 _ground)
    {
        Vector3 currentPos = _testParticle.getPos();
        Vector3 currentVel = _testParticle.getVel();
        float radius = _testParticle.getRadius();
        Vector3 Vel = _testParticle.getVel(), Pos = _testParticle.getPos();
        Vector3 normal = new Vector3(_ground.x,_ground.y,_ground.z);
        normal.Normalize();

        //Test with ground
        if ((currentPos.dot(normal) -_ground.w) < radius)
        {
            //-----------------Calculate Velocity------------------
            //Calculate our new velocity with momentum included
            Vel = -((currentVel.dot(normal))) * normal + (currentVel - (currentVel.cross(normal).cross(normal)));
            Vel+=  m_coefOfRest * currentVel.dot(normal) * normal;

            //If moving parallel to our plane decrease speed due to friction unless it has already stopped
            if(currentVel.sqrMagnitude > 0 && currentVel.dot(normal) == 0)
            {
                Vector3 VelNormalized = Vel;

                VelNormalized.Normalize();
                VelNormalized *=(-1);
                Vector3  friction = new Vector3((1 - normal.x) * VelNormalized.x,(1-normal.y) * VelNormalized.y,(1-normal.z) * VelNormalized.z);
                friction *=((m_coefOfFric*_timeStep)/(_testParticle.getMass()));
                if(friction.length()<Vel.length())
                {
                    Vel +=(friction);
                }
                else if(friction.length()>=Vel.length())
                {
                    Vel.set(0,0,0);
                }
            }

            _testParticle.setVel(Vel);

            //---------------Calculate Position----------------------
            //If particle has a velocity which is not parallel to our plane find its new position
            if(currentVel.length() != 0 || currentVel.dot(normal) != 0)
            {
                Vector3 curPosRadius = currentPos - normal * (radius);
                float t = (_ground.w - normal.dot(curPosRadius)) / normal.dot(normal);
                Vector3 closestPoint = curPosRadius + normal *(t);

                if(t > 0 && t < 1)
                {
                    Pos = closestPoint + normal *(radius);
                }
                else
                {
                    Pos = closestPoint + normal *(radius);
                }
                _testParticle.setPos(Pos);
            }
        }
    }
   
    public static float sqrt(float a)
    {
        return Mathf.Sqrt(a);
    }

    /// @brief Struct to store our wall data neatly
    public struct wall
    {
        public Vector3 center;
        public Vector4 plane;
        public float drawSizeX;
        public float drawSizeY;
        public float drawSizeZ;
        public bool draw;

        public void Center(Vector3 v)
        {
            center = v;
        }

        public void Plane(Vector4 v)
        {
            plane = v;
        }

        public void Draw(bool v)
        {
            draw = v;
        }

        public void DrawSizeX(float v)
        {
            drawSizeX = v;
        }
        public void DrawSizeY(float v)
        {
            drawSizeY = v;
        }
        public void DrawSizeZ(float v)
        {
            drawSizeZ = v;
        }
    }

    private wall Wall;

    private static GameObject m_prefab = null;
    private List<GameObject> m_gos = new List<GameObject>();

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A vector to store all our walls
    public List<wall> m_walls = new List<wall>();

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the coefficient of restitution i.e. how much the particles bouce, it 1 then 0 bounce
    float m_coefOfRest;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the coefficient of friction
    float m_coefOfFric;

}

