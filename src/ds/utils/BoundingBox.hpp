#ifndef BOUNDING_BOX_HPP
#define BOUNDING_BOX_HPP

#ifdef MADVORO_WITH_VCL
    #include <vectorclass.h>
#endif // MADVORO_WITH_VCL

#ifdef MADVORO_WITH_MPI
    #include "mpi/serialize/Serializer.hpp"
#endif // MADVORO_WITH_MPI

#define DIM 3

namespace MadVoro
{
    namespace Geometry
    {
        template<typename T>
        class BoundingBox
                        #ifdef MADVORO_WITH_MPI
                            : public MadVoro::MPI::Serializable
                        #endif // MADVORO_WITH_MPI
        {
        template<typename U>
        friend class BoundingBox;

        private:
            T ll; // leftmost point of the box
            T ur; // rightmost point of the box

            #ifdef MADVORO_WITH_VCL
                Vec4d llVec, urVec;
                Vec4d llPlusUrVec;
                Vec8d boundariesVec;
            #endif // MADVORO_WITH_VCL
            typename T::coord_type width;
            typename T::coord_type widthSquared;

            void recalculateFields()
            {
                typename T::coord_type width = std::max(this->ur[0] - this->ll[0], std::max(this->ur[1] - this->ll[1], this->ur[2] - this->ll[2]));
                this->width = width;
                this->widthSquared = width * width;
                #ifdef MADVORO_WITH_VCL
                    this->llVec = Vec4d(this->ll[0], this->ll[1], this->ll[2], 0);
                    this->urVec = Vec4d(this->ur[0], this->ur[1], this->ur[2], 0);
                    this->llPlusUrVec = this->llVec + this->urVec;
                    this->boundariesVec = Vec8d(this->ll[0], this->ll[1], this->ll[2], this->ur[0], this->ur[1], this->ur[2], 0, 0);
                #endif // MADVORO_WITH_VCL
            }

        public:
            template<typename U>
            inline BoundingBox(const U &_ll, const U &_ur)
            {
                this->ll[0] = _ll[0];
                this->ll[1] = _ll[1];
                this->ll[2] = _ll[2];
                this->ur[0] = _ur[0];
                this->ur[1] = _ur[1];
                this->ur[2] = _ur[2];
                this->recalculateFields();
            };

            inline BoundingBox(): ll(T()), ur(T())
            {
                this->recalculateFields();
            };

            const T &getLL() const{return this->ll;};
            const T &getUR() const{return this->ur;};

            inline void setLL(const T &newLL)
            {
                this->ll = newLL;
                this->recalculateFields();
            }

            inline void setUR(const T &newUR)
            {
                this->ur = newUR;
                this->recalculateFields();
            }
            
            template<typename U>
            inline void setBounds(const U &newLL, const U &newUR)
            {
                this->ll[0] = newLL[0];
                this->ll[1] = newLL[1];
                this->ll[2] = newLL[2];
                this->ur[0] = newUR[0];
                this->ur[1] = newUR[1];
                this->ur[2] = newUR[2];
                this->recalculateFields();
            }

            inline typename T::coord_type getWidth() const
            {
                return this->width;
            }

            inline typename T::coord_type getWidthSquared() const
            {
                return this->widthSquared;
            }

            template<typename U>
            inline bool contains(const U &point) const
            {
                for(int i = 0; i < DIM; i++)
                {
                    if((point[i] < ll[i]) or (point[i] > ur[i]))
                    {
                        return false;
                    }
                }
                return true;
                /*
                Vec8d _point(point[0], point[1], point[2], -point[0], -point[1], -point[2], 0, 0);
                Vec8d _boundaries(this->ll[0], this->ll[1], this->ll[2], -this->ur[0], -this->ur[1], -this->ur[2], 1, 1);
                Vec8db cmp = _point < _boundaries;

                for(int i = 0; i < DIM; i++)
                {
                    if((!cmp[i]) or (!cmp[i + DIM]))
                    {
                        return false;
                    }
                }
                return true;
                */
            }

            /**
             * returns whether the other bounding box is contained in me.
            */
            template<typename U>
            inline bool contained(const BoundingBox<U> &other) const
            {
                // need: other.ll[i] >= this->ll[i] and other.ur[i] <= this->ur[i] for all i
                for(int i = 0; i < DIM; i++)
                {
                    if((other.ll[i] < this->ll[i]) or (other.ur[i] > this->ur[i]))
                    {
                        return false;
                    }
                }
                return true;
            }

            template<typename U>
            inline bool intersects(const BoundingBox<U> &other) const
            {
                return this->contains(other.getLL()) or this->contains(other.getUR()) or other.contains(this->ll) or other.contains(this->ur);
            }

            template<typename U>
            T closestPoint(const U &point) const;

            template<typename U>
            T closestPointToOther(const BoundingBox<U> &otherBox) const;
            
            template<typename U>
            inline typename T::coord_type distanceSquared(const U &point) const
            {
                T closestPoint = this->closestPoint(point);
                #ifdef MADVORO_WITH_VCL
                    Vec4d closestPointVec(closestPoint[0] - point[0], closestPoint[1] - point[1], closestPoint[2] - point[2], 0);
                    Vec4d pointSquaredVec = closestPointVec * closestPointVec;
                    return pointSquaredVec[0] + pointSquaredVec[1] + pointSquaredVec[2];
                #else
                    typename T::coord_type closestDistance = 0;
                    for(int i = 0; i < DIM; i++)
                    {
                        closestDistance += (point[i] - closestPoint[i]);
                    }
                    return closestDistance;
                #endif // MADVORO_WITH_VCL
            }

            template<typename U>
            inline T furthestPoint(const U &point) const
            {
                #ifdef MADVORO_WITH_VCL
                    Vec4d _point(point[0], point[1], point[2], 0);
                    Vec4db cmp = (2 * _point) > this->llPlusUrVec;
                #endif // MADVORO_WITH_VCL

                T furthestPoint;

                for(int i = 0; i < DIM; i++)
                {
                    #ifdef MADVORO_WITH_VCL
                        if(cmp[i])
                    #else // MADVORO_WITH_VCL
                        if(point[i] < ((this->ll[i] + this->ur[i]) * 0.5))
                    #endif // MADVORO_WITH_VCL
                    {
                        // that means point[i] > (this->ll[i] + this->ur[i]) / 2
                        furthestPoint[i] = this->ll[i];
                    }
                    else
                    {
                        furthestPoint[i] = this->ur[i];
                    }
                }
                return furthestPoint;
            }

            friend inline std::ostream &operator<<(std::ostream &stream, const BoundingBox<T> &box)
            {
                return stream << "BoundingBox(" << box.ll << ", " << box.ur << ")";
            }

            #ifdef MADVORO_WITH_MPI

            force_inline size_t dump(MadVoro::MPI::Serializer *serializer) const override
            {
                size_t bytes = 0;
                bytes += this->ll.dump(serializer);
                bytes += this->ur.dump(serializer);
                return bytes;
            }
            
            force_inline size_t load(const MadVoro::MPI::Serializer *serializer, size_t byteOffset) override
            {
                size_t bytesRead = 0;
                bytesRead += this->ll.load(serializer, byteOffset);
                bytesRead += this->ur.load(serializer, byteOffset + bytesRead);
                this->recalculateFields();
                return bytesRead;
            }

            #endif // MADVORO_WITH_MPI

            template<typename U>
            inline bool operator==(const BoundingBox<U> &other)
            {
                return ((this->ll == other.ll) and (this->ur == other.ur));
            }
        };

        template<typename T>
        template<typename U>
        T BoundingBox<T>::closestPoint(const U &point) const
        {
            #ifdef MADVORO_WITH_VCL
                Vec8d _point(point[0], point[1], point[2], point[0], point[1], point[2], 0, 0);
                Vec8db cmp = _point < this->boundariesVec; // _point < _boundaries;
            #endif // MADVORO_WITH_VCL
            T closestPoint;
            for(int i = 0; i < DIM; i++)
            {
                #ifdef MADVORO_WITH_VCL
                    if(cmp[i])
                #else
                    if(point[i] < this->ll[i])
                #endif // MADVORO_WITH_VCL
                {
                    // that means point[i] < this->ll[i]
                    closestPoint[i] = this->ll[i];
                }
                else
                {
                    #ifdef MADVORO_WITH_VCL
                        if(cmp[i + DIM])
                    #else // MADVORO_WITH_VCL
                        if(point[i] < this->ur[i])
                    #endif // MADVORO_WITH_VCL
                    {
                        // that means point[i] < this->ur[i]
                        closestPoint[i] = point[i];
                    }
                    else
                    {
                        closestPoint[i] = this->ur[i];
                    }
                }
            }
            return closestPoint;
        }

        template<typename T>
        template<typename U>
        T BoundingBox<T>::closestPointToOther(const BoundingBox<U> &otherBox) const
        {
            // TODO: vectorization!
            T closestPoint;
            for(int i = 0; i < DIM; i++)
            {
                // check where's ll[i] and ur[i] compared to ll2[i] and ur2[i]
                typename T::coord_type ll1 = ll[i], ur1 = ur[i], ll2 = otherBox.ll[i], ur2 = otherBox.ur[i];
                if(ll2 > ur1)
                {
                    closestPoint[i] = ur1;
                }
                else if(ur2 < ll1)
                {
                    closestPoint[i] = ll1;
                }
                else
                {
                    // ll1 < ll2 < ur1
                    closestPoint[i] = ll2;
                }
            }
            return closestPoint;
        }
    }
}

#endif // BOUNDING_BOX_HPP