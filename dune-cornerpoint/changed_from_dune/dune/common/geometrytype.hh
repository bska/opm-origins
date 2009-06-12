#ifndef DUNE_GEOMETRY_TYPE_HH
#define DUNE_GEOMETRY_TYPE_HH

/** \file
    \brief A unique label for each type of element that can occur in a grid
*/

#include <dune/common/exceptions.hh>

namespace Dune {

    /** \brief Unique label for each type of entities that can occur in DUNE grids

    This class has to be extended if a grid implementation with new entity types
    is added to DUNE.
	
	\ingroup COGeometryType	
     */
    class GeometryType
    {
    public:
        /** \brief Each entity can be tagged by one of these basic types 
            plus its space dimension */
        enum BasicType {
		  simplex,   //!< Simplicial element in any nonnegative dimension
		  cube,      //!< Cube element in any nonnegative dimension
		  pyramid,   //!< Four sided pyramid in three dimensions
		  prism,     //!< Prism element in three dimensions
		  singular   //!< Singular element, i.e. element with only one point represented (centroid)
		};

    private:

        /** \brief Basic type of the element */
        BasicType basicType_;

        /** \brief Dimension of the element */
        unsigned int dim_;

    public:
        /** \brief Default constructor, not initializing anything */
        GeometryType ()
 {}

        /** \brief Constructor */
        GeometryType(BasicType basicType, unsigned int dim)
            : basicType_(basicType), dim_(dim)
        {}

        /** \brief Constructor for vertices and segments
            \todo Add check for dim={0,1} when compiled with a suitable flag
        */
        explicit GeometryType(unsigned int dim)
            : basicType_(cube), dim_(dim)
        {}

        /** @name Setup Methods */
        /*@{*/

        /** \brief Make a vertex */
        void makeVertex() {dim_ = 0;}
        
        /** \brief Make a line segment */
        void makeLine() {dim_ = 1;}
        
        /** \brief Make a triangle */
        void makeTriangle() {basicType_ = simplex; dim_ = 2;}
        
        /** \brief Make a quadrilateral */
        void makeQuadrilateral() {basicType_ = cube; dim_ = 2;}
        
        /** \brief Make a tetrahedron */
        void makeTetrahedron() {basicType_ = simplex; dim_ = 3;}
        
        /** \brief Make a pyramid */
        void makePyramid() {basicType_ = pyramid;}
        
        /** \brief Make a prism */
        void makePrism() {basicType_ = prism;}
        
        /** \brief Make a hexahedron */
        void makeHexahedron() {basicType_ = cube; dim_ = 3;}
        
        /** \brief Make a simplex of given dimension */
        void makeSimplex(unsigned int dim) {basicType_ = simplex; dim_ = dim;}
        
        /** \brief Make a hypercube of given dimension */
        void makeCube(unsigned int dim) {basicType_ = cube; dim_ = dim;}
        
        /** \brief Make a singular of given dimension */
        void makeSingular(unsigned int dim) {basicType_ = singular; dim_ = dim;}
        
        /*@}*/

        
        /** @name Query Methods */
        /*@{*/
        /** \brief Return true if entity is a vertex */
        bool isVertex() const {return dim_==0;}
        
        /** \brief Return true if entity is a line segment */
        bool isLine() const {return dim_==1;}
        
        /** \brief Return true if entity is a triangle */
        bool isTriangle() const {return basicType_==simplex && dim_==2;}
        
        /** \brief Return true if entity is a quadrilateral */
        bool isQuadrilateral() const {return basicType_==cube && dim_==2;}
        
        /** \brief Return true if entity is a tetrahedron */
        bool isTetrahedron() const {return basicType_==simplex && dim_==3;}
        
        /** \brief Return true if entity is a pyramid */
        bool isPyramid() const {return basicType_==pyramid;}
        
        /** \brief Return true if entity is a prism */
        bool isPrism() const {return basicType_==prism;}
        
        /** \brief Return true if entity is a hexahedron */
        bool isHexahedron() const {return basicType_==cube && dim_==3;}
        
        /** \brief Return true if entity is a simplex of any dimension */
        bool isSimplex() const {return basicType_==simplex || dim_ < 2;}
        
        /** \brief Return true if entity is a cube of any dimension */
        bool isCube() const {return basicType_==cube || dim_ < 2;}
        
        /** \brief Return true if entity is a singular of any dimension */
        bool isSingular() const {return basicType_==singular;}
        
        /** \brief Return dimension of the entity */
        unsigned int dim() const {return dim_;}

        /** \brief Return the basic type of the entity */
        BasicType basicType() const {return basicType_;}

        /*@}*/

        /** \brief Check for equality. This method knows that in dimension 0 and 1
			all BasicTypes are equal.
		 */
        bool operator==(const GeometryType& other) const {
            return ( (dim()==0 && other.dim()==0) 
                     || (dim()==1 && other.dim()==1)
                     || (dim()==other.dim() && basicType_==other.basicType_) );
        }

        /** \brief Check for inequality */
        bool operator!=(const GeometryType& other) const {
            return ! ((*this)==other);
        }

        /** \brief lesser operation for use with maps */
        bool operator < (const GeometryType& other) const {	    
            if (dim() != other.dim())
                return dim() < other.dim();
            else if (dim()==0 || dim()==1)
                return false;
            
            return basicType_ < other.basicType_;
        }

        /** \brief Prints the type to an output stream */
        friend std::ostream& operator<< (std::ostream& s, const GeometryType& a)
        {
            switch (a.basicType_) {
            case simplex:   
                s << "(simplex, " << a.dim_ << ")";
                break;
            case cube:      
                s << "(cube, " << a.dim_ << ")";
                break;
            case pyramid:
                s << "pyramid"; 
                break;
            case prism:
                s << "prism";
		break;
	    case singular:
		s << "singular";
		break;
            }
            
            return s;
        }

    };

    /** \brief Prints a GeometryType::BasicType to an output stream */
    inline std::ostream& operator<< (std::ostream& s, GeometryType::BasicType type)
    {
      switch (type) {
      case GeometryType::simplex: s << "simplex"; break;
      case GeometryType::cube:    s << "cube";    break;
      case GeometryType::pyramid: s << "pyramid"; break;
      case GeometryType::prism:   s << "prism";   break;
      default: s << "[unknown GeometryType::BasicType: " << int(type) << "]";
      }
      return s;
    }
}

#endif
