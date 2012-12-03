// -*- mode: c++; mode: rebox; fill-column: 100; rebox-min-fill-column: 100 -*-

/******************************************************************************
 * src/external/tuples.h
 *
 * Utility header to implement tuples and further tools for skew3.
 *
 * Cleaned up and parts reimplemented by Timo Bingmann (2012)
 * Originally implemented in STXXL by Jens Mehnert (2004)
 *
 ******************************************************************************
 * Copyright (C) 2004 Jens Mehnert
 * Copyright (C) 2012 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

namespace tuples {

/**
 * Comparator class for ordering of pairs by first component.
 */
template <typename Type>
struct less
{
    typedef Type       value_type;

    bool operator()(const value_type& a, const value_type& b) const {
        return (a < b);
    }

    static value_type min_value() { return std::numeric_limits<value_type>::min(); }
    static value_type max_value() { return std::numeric_limits<value_type>::max(); }
};

/**
 * A pair class, that is ostream-able!
 * @param T1 first component type.
 * @param T2 second component type.
 */

template<class T1, class T2>
class pair
{
public:
    typedef T1	first_type;
    typedef T2	second_type;

    T1		first;
    T2		second;

    pair()
        : first(T1()), second(T2())
    {}

    pair(const T1& a,const T2& b)
        : first(a), second(b)
    {}

    pair(const pair<T1,T2>& t)
        : first(t.first), second(t.second)
    {}

    pair& operator=(const pair<T1,T2>& t)
    {
        first = t.first; second = t.second;
        return *this;
    }

    bool operator==(const pair<T1,T2>& x)
    {
        return (x.first == first) && (x.second == second);
    }

    bool operator!=(const pair<T1,T2>& x)
    {
        return (x.first != first) && (x.second != second);
    }

    static std::ostream& description (std::ostream& os)
    {
        return os << "(first,second)";
    }

    friend std::ostream& operator<< (std::ostream& os, const pair& p)
    {
        return os << "(" << p.first  << "," << p.second << ")";
    }

    static pair min_value() { return pair(std::numeric_limits<first_type>::min(),
                                          std::numeric_limits<second_type>::min()); }

    static pair max_value() { return pair(std::numeric_limits<first_type>::max(),
                                          std::numeric_limits<second_type>::max()); }
};

/**
 * Comparator class for ordering of pairs by first component.
 */
template <typename PairType>
struct pair_less1st
{
    typedef PairType       value_type;

    bool operator()(const value_type& a, const value_type& b) const {
        if (a.first == b.first)
            return (a.second < b.second);
        return (a.first < b.first);
    }

    static value_type min_value() { return value_type::min_value(); }
    static value_type max_value() { return value_type::max_value(); }
};

/**
 * Comparator class for ordering of pairs by first component.
 */
template <typename PairType>
struct pair_greater1st
{
    typedef PairType       value_type;

    bool operator()(const value_type& a, const value_type& b) const {
        if (a.first == b.first)
            return (a.second > b.second);
        return (a.first > b.first);
    }

    static value_type min_value() { return value_type::max_value(); }
    static value_type max_value() { return value_type::min_value(); }
};

/**
 * Comparator class for ordering of pairs by second component.
 */
template <typename PairType>
struct pair_less2nd
{
    typedef PairType       value_type;

    bool operator()(const value_type& a, const value_type& b) const {
        if (a.second == b.second)
            return (a.first < b.first);
        return (a.second < b.second);
    }

    static value_type min_value() { return value_type::min_value(); }
    static value_type max_value() { return value_type::max_value(); }
};

/**
 * A triple class.
 * @param T1 first component type.
 * @param T2 second component type.
 * @param T3 third component type.
 */

template<class T1, class T2, class T3>
class triple
{
public:
    typedef T1	first_type;
    typedef T2	second_type;
    typedef T3	third_type;

    T1		first;
    T2		second;
    T3		third;

    triple()
        : first(T1()), second(T2()), third(T3())
    {}

    triple(const T1& a,const T2& b,const T3& c)
        : first(a), second(b), third(c)
    {}

    triple(const triple<T1,T2,T3>& t)
        : first(t.first), second(t.second), third(t.third)
    {}

    triple& operator=(const triple<T1,T2,T3>& t)
    {
        first = t.first; second = t.second; third = t.third;
        return *this;
    }

    bool operator==(const triple<T1, T2, T3>& x)
    {
        return (x.first == first) && (x.second == second) && (x.third == third);
    }

    bool operator!=(const triple<T1, T2, T3>& x)
    {
        return (x.first != first) && (x.second != second) && (x.third != third);
    }

    static std::ostream& description (std::ostream& os)
    {
        return os << "(first,second,third)";
    }

    friend std::ostream& operator<< (std::ostream& os, const triple& t)
    {
        return os << "(" << t.first  << "," << t.second << "," << t.third << ")";
    }

    static triple min_value() { return triple(std::numeric_limits<first_type>::min(),
                                              std::numeric_limits<second_type>::min(),
                                              std::numeric_limits<third_type>::min()); }

    static triple max_value() { return triple(std::numeric_limits<first_type>::max(),
                                              std::numeric_limits<second_type>::max(),
                                              std::numeric_limits<third_type>::max()); }
};

template <class Triple>
struct triple_less1st
{
    typedef Triple	value_type;

    bool operator() (const value_type& a, const value_type& b) const {
        return a.first < b.first;
    }

    static value_type min_value() { return value_type::min_value(); }
    static value_type max_value() { return value_type::max_value(); }
};

/**
 * A quadruple class.
 * @param T1 first component type.
 * @param T2 second component type.
 * @param T3 third component type.
 * @param T4 fourth component type.
 */

template<class T1, class T2, class T3, class T4>
class quad
{
public:
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;
    typedef T4 fourth_type;

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
    //bool fin;

    quad()
        : first(T1()), second(T2()), third(T3()), fourth(T4())
    {}

    quad(const T1& a, const T2& b, const T3& c, const T4& d)
        : first(a), second(b), third(c), fourth(d)
    {
    }

    quad(const quad<T1,T2,T3,T4>& t)
        : first(t.first), second(t.second), third(t.third), fourth(t.fourth)
    {}

    quad& operator=(const quad<T1,T2,T3,T4>& t)
    {
        first=t.first; second=t.second; third=t.third; fourth=t.fourth;
        return *this;
    }

    bool operator==(const quad<T1,T2,T3,T4>& x)
    {
        return (x.first == first) && (x.second == second) && (x.third == third) && (x.fourth== fourth);
    }

    bool operator!=(const quad<T1,T2,T3,T4>& x)
    {
        return (x.first != first) && (x.second != second) && (x.third != third) && (x.fourth != fourth);
    }

    const T1& get_1() const { return first; }
    const T2& get_2() const { return second; }
    const T3& get_3() const { return third; }
    const T4& get_4() const { return fourth; }

    void set_1(const T1& x) { first = x; }
    void set_2(const T2& x) { second = x; }
    void set_3(const T3& x) { third = x; }
    void set_4(const T4& x) { fourth = x; }

    friend std::ostream& operator<< (std::ostream& os, const quad& q)
    {
        return os << "(" << q.first  << "," << q.second << "," << q.third << "," << q.fourth << ")";
    }

    static quad min_value() { return quad(std::numeric_limits<first_type>::min(),
                                          std::numeric_limits<second_type>::min(),
                                          std::numeric_limits<third_type>::min(),
                                          std::numeric_limits<fourth_type>::min()); }

    static quad max_value() { return quad(std::numeric_limits<first_type>::max(),
                                          std::numeric_limits<second_type>::max(),
                                          std::numeric_limits<third_type>::max(),
                                          std::numeric_limits<fourth_type>::max()); }
};

/**
 * A quintuple class.
 * @param T1 first component type.
 * @param T2 second component type.
 * @param T3 third component type.
 * @param T4 fourth component type.
 * @param T5 fifth component type.
 */
template<class T1, class T2, class T3, class T4,class T5>
class quint
{
public:
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;
    typedef T4 fourth_type;
    typedef T5 fifth_type;

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
    T4 fifth;
    //bool fin;

    quint()
        : first(T1()), second(T2()), third(T3()), fourth(T4()), fifth(T5())
    {}

    quint(const T1& a,const T2& b,const T3& c, const T4& d,const T5& e)
        : first(a), second(b), third(c), fourth(d), fifth(e)
    {}

    quint(const quint<T1,T2,T3,T4,T5>& t)
        : first(t.first), second(t.second), third(t.third), fourth(t.fourth), fifth(t.fifth)
    {}

    quint& operator=(const quint<T1,T2,T3,T4,T5>& t)
    {
        first=t.first; second=t.second; third=t.third; fourth=t.fourth; fifth=t.fifth;
        return *this;
    }

    bool operator==(const quint<T1,T2,T3,T4,T5>& x)
    {
        return x.first == first && x.second == second && x.third == third && x.fourth== fourth && x.fifth== fifth;
    }

    bool operator!=(const quint<T1,T2,T3,T4,T5>& x){
        return x.first != first && x.second != second && x.third != third && x.fourth != fourth && x.fifth != fifth;
    }

    const T1& get_1() const { return first; }
    const T2& get_2() const { return second; }
    const T3& get_3() const { return third; }
    const T4& get_4() const { return fourth; }
    const T5& get_5() const { return fifth; }

    void set_1(const T1& x) { first = x; }
    void set_2(const T2& x) { second = x; }
    void set_3(const T3& x) { third = x; }
    void set_4(const T4& x) { fourth = x; }
    void set_5(const T5& x) { fifth = x; }

    friend std::ostream& operator<< (std::ostream& os, const quint& q)
    {
        return os << "(" << q.first  << "," << q.second << "," << q.third << "," << q.fourth << "," << q.fifth << ")";
    }

    static quint min_value() { return quint(std::numeric_limits<first_type>::min(),
                                            std::numeric_limits<second_type>::min(),
                                            std::numeric_limits<third_type>::min(),
                                            std::numeric_limits<fourth_type>::min(),
                                            std::numeric_limits<fifth_type>::min()); }

    static quint max_value() { return quint(std::numeric_limits<first_type>::max(),
                                            std::numeric_limits<second_type>::max(),
                                            std::numeric_limits<third_type>::max(),
                                            std::numeric_limits<fourth_type>::max(),
                                            std::numeric_limits<fifth_type>::max()); }
};

/* Iteratoradapters from Roman */

template <class InputIterator_>
class it_rg
{
    InputIterator_ current_, end_;
    it_rg(){}
public:

    typedef typename std::iterator_traits<InputIterator_>::value_type value_type;

    it_rg(InputIterator_ begin, InputIterator_ end):current_(begin),end_(end){}

    it_rg(const it_rg& a):current_(a.current_),end_(a.end_){}

    const value_type& operator *()const{
        assert(end_ !=current_);
        return *current_;
    }

    it_rg& operator ++(){
        assert(end_ != current_);
        ++current_;
        return *this;
    }

    bool empty() const{
        return (current_ == end_ );
    }
};

template<class InputIterator>
it_rg<InputIterator> create_ls(InputIterator begin, InputIterator end){
  return it_rg<InputIterator>(begin,end);
}

template<class OutputIterator_ , class LazyAlgorithm_>
OutputIterator_ flush_2_it(OutputIterator_ out, LazyAlgorithm_ & in)
{
    while(!in.empty()){
        *out = *in;
        ++out;
        ++in;
    }
    return out;
}

/* Output the inputstream for debug issues */

template<class A1>
class s_print
{
public:
    typedef typename A1::value_type value_type;

private:
    A1& A;
    value_type curr;
    unsigned int counter;
    bool p_flag;

public:
    s_print(A1& x, bool p)
        : A(x),
          counter(0),
          p_flag(p)
    {
        assert(!A.empty());
        curr = *A;
    }

    const value_type& operator*() const
    {
        if (p_flag) std::cout << "[" << counter << "] = " << curr << std::endl;
        return curr;
    }

    s_print& operator++()
    {
        ++A; ++counter;
        if (!A.empty()) curr = *A;
        return *this;
    }

    bool empty() const
    {
        if (A.empty() && p_flag) std::cout << std::endl;
        return A.empty();
    }
};

/**
 * Counter for creating tuple indexes for example.
 */
template <class ValueType>
struct counter
{
    typedef ValueType value_type;

    value_type cnt;

    counter() : cnt(0) {}

    const value_type & operator*() const
    {
        return cnt;
    }

    counter& operator++ ()
    {
        ++cnt;
        return *this;
    }

    bool empty() const
    {
        return false;
    }
};

/**
 * Wrap a stream with a debug output layer
 */
template <bool dbgflag, typename Parent>
class print : public Parent
{
private:
    counter<size_t> cnt;

public:

    template <typename A>
    print(A& a) : Parent(a) {}

    template <typename A, typename B>
    print(A& a, B& b) : Parent(a,b) {}

    template <typename A, typename B, typename C>
    print(A& a, B& b, C& c) : Parent(a,b,c) {}

    template <typename A, typename B, typename C, typename D>
    print(A& a, B& b, C& c, D& d) : Parent(a,b,c,d) {}

    print& operator++()
    {
        if (*cnt == 0 && !Parent::empty())
            DBG(dbgflag, *cnt << " = " << Parent::operator*());

        Parent::operator++();
        ++cnt;

        if (!Parent::empty())
            DBG(dbgflag, *cnt << " = " << Parent::operator*());

        return *this;
    }
};

/**
 * A four-tuple class for doubling.
 * @param T1 first component type.
 * @param T2 second/third fourth/fifth component type.
 */

template<class T1, class T2>
class doubling_ntuple
{
public:
    typedef T1 first_type;
    typedef T2 second_type;

    T1 first;
    T2 second;
    T2 third;
    T2 fourth;
    T2 fifth;

    doubling_ntuple()
        : first(T1()), second(T2()), third(T2()), fourth(T2()), fifth(T2())
    {}

    doubling_ntuple(const T1& a, const T2& c, const T2& d, const T2& e, const T2& f)
        : first(a), second(c), third(d), fourth(e), fifth(f)
    {}

    doubling_ntuple(const doubling_ntuple<T1,T2>& t)
        : first(t.first), second(t.second), third(t.third), fourth(t.fourth), fifth(t.fifth)
    {}

    doubling_ntuple& operator=(const doubling_ntuple<T1,T2>& t)
    {
        first = t.first; second = t.second; third = t.third; fourth = t.fourth; fifth = t.fifth;
        return *this;
    }

    bool operator==(const doubling_ntuple<T1,T2>& x)
    {
        return (x.first == first) && (x.second == second) && (x.third == third) && (x.fourth == fourth) && (x.fifth == fifth);
    }

    bool operator!=(const doubling_ntuple<T1, T2>& x)
    {
        return (x.first != first) && (x.second != second) && (x.third != third) && (x.fourth != fourth) && (x.fifth != fifth);
    }

    /* To encode the finish bit in the second component, i need its absolute value */

    /* Getter */

    const T1& value_1() const { return first; }
    const T2& value_2() const { return second; }
    const T2& value_3() const { return third; }
    const T2& value_4() const { return fourth; }
    const T2& value_5() const { return fifth; }

    /* Setter */

    void set_1(const T1& x) { first = x; }
    void set_2(const T2& x) { second = x; }
    void set_3(const T2& x) { third = x; }
    void set_4(const T2& x) { fourth = x; }
    void set_5(const T2& x) { fifth = x; }

    void to_cout() const {
        std::cout << "(" << first  << "," << second << "," << third << "," << fourth<< "," << fifth << ") " << std::endl;
    }
};

} // namespace tuples
