namespace octet {
  namespace scene {
    class mesh_sprite : public mesh {
      mat4t transform;
      quad _quad;

    public:
      RESOURCE_META (mesh_sprite)
        mesh_sprite () {}

      /*  todo: override: 
      *   virtual btCollisionShape *get_bullet_shape() - mesh member   
      */ 
    
      mesh_sprite (vec3 pos, vec2 size, mat4t transform) {
        //todo: remove normals
        set_default_attributes ();
        this->transform = transform;
        _quad = quad (pos, size / 2);
        vec3 corners[4] = {
          vec3(pos+size/2),
          vec3(pos-size/2),
          vec3(pos+size.rot90()/2),
          vec3(pos-size.rot90()/2)
        };
        set_aabb (aabb (corners, corners+3));
        update ();
      }

      /// Generate mesh from parameters.
      virtual void update () {
        aabb aabb_ = get_aabb ();
        mesh::set_shape<math::aabb, mesh::vertex> (aabb_, transform, 1);
        //dump(log("zzz\n"));
      }

      /// Serialise the box
      void visit (visitor &v) {
        mesh::visit (v);
        if ( v.is_reader () ) {
          update ();
        }
      }

/*
*     todo: add add this later when preparing game physics. 
*
#ifdef OCTET_BULLET
      /// Get a bullet shape object for this mesh
      btCollisionShape *get_bullet_shape () {
        return new btBoxShape (get_btVector3 (get_aabb ().get_half_extent ()));
      }

      /// Get a static bullet shape object for this mesh
      btCollisionShape *get_static_bullet_shape () {
        return new btBoxShape (get_btVector3 (get_aabb ().get_half_extent ()));
      }
#endif
*/
    };
  }
}

