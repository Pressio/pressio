
#ifndef COMMONBASE_SAME_INSTANCE_HPP_INCLUDED
#define COMMONBASE_SAME_INSTANCE_HPP_INCLUDED


namespace commonBase{
namespace details{
namespace meta {

template<class T1 , class T2>
struct same_instance_impl
{ 
    static bool same_instance( const T1& /* x1 */ , const T2& /* x2 */ )
    {
        return false;
    }
};

template< class T >
struct same_instance_impl<T,T>
{ 
    static bool same_instance( const T &x1 , const T &x2 )
    {
        return (&x1 == &x2);
    }
};


template< class T1 , class T2 >
bool same_instance( const T1 &x1 , const T2 &x2 )
{
    return same_instance_impl< T1 , T2 >::same_instance( x1 , x2 );
}


} // namespace meta
} // namespace details
} // namespace commonBase

#endif
