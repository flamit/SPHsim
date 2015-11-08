using UnityEngine;
//----------------------------------------------------------------------------------------------------------------------
/// @file particle.h
/// particle class, it will store all the data for my particles for my simulation

//----------------------------------------------------------------------------------------------------------------------


//#include <ngl/Vec4.h>
//#include <ngl/Camera.h>
//#include <ngl/TransformStack.h>
//#include <ngl/VAOPrimitives.h>
//#include <ngl/ShaderLib.h>


public class particle
{

    /// @brief our constructor to create our particle
    /// @param _pos -  a Vector3* of the position of our particle
    /// @param _mass -  the mass of our particle used for later calculations
    /// @param _partOfMesh - if part of mesh will store a pointer to the vector location else will create a local one
    public particle(Vector3 _pos, float _mass)
    {
        m_position = _pos;
        m_mass = _mass;
    }

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A function to return our particles position
    public Vector3 getPos()  {return m_position;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A function to return our particles Acceleration
    public Vector3 getAcc() {return m_acceleration;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A function to return our particles velocity
    public Vector3 getVel() {return m_velocity;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A function to return the radius of our particle
    public float getRadius() {return m_radius;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A function to update our Postion once calculated
    /// @param _pos Vertex position
    public void setPos(Vector3 _pos)
    {
        if (float.IsNaN(_pos.x) || float.IsNaN(_pos.y) || float.IsNaN(_pos.z))
        {
            return;
        }

        m_position = _pos;
    }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A function to update our Acceleration once calculated
    /// @param _acc the acceleration Vertext
    public   void setAcc(Vector3 _acc){m_acceleration = _acc;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A function to update our Velocity once calculated
    /// @param _vel the Velocity Vertex
    public   void setVel(Vector3 _vel){m_velocity = _vel;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the mass of our particle
    /// @param _mass - the mass
    public   void setMass(float _mass){m_mass = _mass;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief returns the mass of our particle
    public   float getMass(){return m_mass;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the dynamic density of our particle
    /// @param _density - the density
    public   void setDensity(float _density){m_density = _density;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief returns the rest density of our particle
    public   float getRestDensity(){return m_restDensity;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the rest density of our particle
    /// @param _density - the rest density
    public   void setRestDensity(float _density){m_restDensity = _density;}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief returns the density of our particle
    public   float getDensity(){return m_density;}

    //helper vars to draw primitive circles in unity opengl

    public static int lineCount = 3;
    public static Material lineMaterial;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief draws our particles
    /// @param trans - our scenes transform stack to get our M matrix
    public void draw(Camera cam, Transform tr)
    {
        if (float.IsInfinity(m_position.x) || float.IsInfinity(m_position.y) || float.IsInfinity(m_position.z))
        {
            return;
        }
        tr.position = m_position;
        lineMaterial.SetPass(0);
        GL.PushMatrix();
        // Set transformation matrix for drawing to
        // match our transform
        GL.MultMatrix(tr.localToWorldMatrix);
        GL.Color(new Color(0, (255f / m_density) * 100f, (255f / m_density) * 100f));
        // Draw lines
        GL.Begin(GL.LINES);
        for (int i = 0; i < lineCount; ++i)
        {
            float a = i / (float) lineCount;
            float prevAngle = (i - 1) / (float) lineCount * (Mathf.PI * 2f);
            float angle = a * Mathf.PI * 2;
            GL.Vertex3(Mathf.Cos(prevAngle) * m_radius, Mathf.Sin(prevAngle) * m_radius, 0);
            GL.Vertex3(Mathf.Cos(angle) * m_radius, Mathf.Sin(angle) * m_radius, 0);
        }
        GL.End();
        GL.PopMatrix();
    }

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A variable to store our particles position
    /// @brief this is a pointer so that we can update the postion of the mesh data when we do our calculations rather then having to copy it over every time
    private Vector3 m_position =  Vector3.zero;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A variable to store our particles mass
    private float m_mass;

    /// @brief A variable to store our particles acceleration
    private Vector3 m_acceleration =  Vector3.zero;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A variable to store our particles velocity
    private Vector3 m_velocity =  Vector3.zero;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief A variable to store the dynamic density of our particle
    private float m_density;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief The resting density of our particle
    private float m_restDensity;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief The radius of our particle for collision detection
    private static float m_radius = 0.1f;
}
