 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: index_pizzachili.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_PIZZACHILI_H
#define SEQAN_HEADER_INDEX_PIZZACHILI_H

#include <seqan/index/pizzachili_api.h>
#include <seqan/index/index_pizzachili_string.h>

namespace SEQAN_NAMESPACE_MAIN {

/**
.Tag.Pizza & Chili Index Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Spec.Pizza & Chili Index@ index.
..cat:Index
..tag.PizzaChili_Text:The original text the index is based on.
..tag.PizzaChili_Compressed:The compressed suffix array.
...remarks:Pizza & Chili indices are compressed indices. Hence, this fibre is used for searching in the index.
..see:Metafunction.Fibre
..see:Function.getFibre
..see:Spec.Pizza & Chili Index
*/

struct _Fibre_PizzaChili_Compressed;

typedef Tag<_Fibre_Text> const Fibre_PizzaChili_Text;
typedef Tag<_Fibre_PizzaChili_Compressed> const Fibre_PizzaChili_Compressed;

typedef Fibre_PizzaChili_Text PizzaChili_Text;
typedef Fibre_PizzaChili_Compressed PizzaChili_Compressed;

/**
.Spec.Pizza & Chili Index:
..summary:An adapter for the Pizza & Chili index API.
..remarks:
..cat:Index
..general:Class.Index
..signature:Index<TText, PizzaChili<TSpec> >
..param.TText:The text type.
...type:Class.String
..param.TSpec:Tag specifying the Pizza & Chili index library to use.
...type.Tag.Pizza & Chili Index Tags
..see:Spec.Pizza & Chili String
..see:Tag.Pizza & Chili Index Fibres
..see:Tag.Index Find Algorithm
*/

template <typename TText, typename TSpec>
class Index<TText, PizzaChili<TSpec> > {
public:
    typedef typename Value<TText>::Type TValue;
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;

    impl::index_t index_handle;
    Holder<String<TValue, PizzaChili<TSpec> > > text;

    Index() : index_handle(0), text() { }

    Index(Index& other) : index_handle(0), text() {
SEQAN_CHECKPOINT
        // Explicitly request the other's index text.
        setIndexText(*this, indexText(other));
    }

    Index(Index const& other) : index_handle(0), text() {
SEQAN_CHECKPOINT
        // Explicitly request the other's index text.
        setIndexText(*this, indexText(other));
    }

    template <typename TOtherText>
    Index(TOtherText& txt) : index_handle(0), text() {
SEQAN_CHECKPOINT
        setIndexText(*this, txt);
    }

    ~Index() {
SEQAN_CHECKPOINT
        clear(*this);
    }

    Index& operator =(Index const& other) {
SEQAN_CHECKPOINT
        if (this == &other)
            return *this;

        clear(*this);
        setIndexText(*this, indexText(other));

        return *this;
    }

private:
    //Index(Index const& other);
    //Index operator =(Index const& other);
};

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
struct Fibre<Index<TText, PizzaChili<TSpec> >, PizzaChili_Text> {
    typedef String<typename Value<TText>::Type, PizzaChili<TSpec> >& Type;
};

template <typename TText, typename TSpec>
struct Fibre<Index<TText, PizzaChili<TSpec> > const, PizzaChili_Text> {
    typedef String<typename Value<TText>::Type, PizzaChili<TSpec> > const& Type;
};

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TText, typename TSpec>
    inline void
    clearIndex(Index<TText, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
        typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;

        if (me.index_handle != 0) {
            impl::error_t e =
                TCodeProvider::free_index(me.index_handle);

            if (e != 0)
                SEQAN_ABORT(TCodeProvider::error_index(e));

            me.index_handle = 0;
        }
    }
} // namespace impl

template <typename TText, typename TSpec>
inline void
clear(Index<TText, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
    impl::clearIndex(me);
    clear(me.text);
}

//////////////////////////////////////////////////////////////////////////////

namespace impl {

    template <typename TText, typename TSpec>
    inline char const*
    getOptionsString(Index<TText, PizzaChili<TSpec> >& /*me*/) {
SEQAN_CHECKPOINT
        return "";
    }

    template <typename TText>
    inline char const*
    getOptionsString(Index<TText, PizzaChili<PizzaChili_SA> >& /*me*/) {
SEQAN_CHECKPOINT
        return "copy_text";
    }

    template <typename TText>
    inline char const*
    getOptionsString(Index<TText, PizzaChili<PizzaChili_FM> >& /*me*/) {
SEQAN_CHECKPOINT
        return "-a 0";
    }

    template <typename TText>
    inline char const*
    getOptionsString(Index<TText, PizzaChili<PizzaChili_RSA> >& /*me*/) {
SEQAN_CHECKPOINT
        return "copy_text";
    }

} // namespace impl

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> >, PizzaChili_Text>::Type
indexText(Index<TText, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
    return getFibre(me, PizzaChili_Text());
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> > const, PizzaChili_Text>::Type
indexText(Index<TText, PizzaChili<TSpec> > const& me) {
SEQAN_CHECKPOINT
    return getFibre(me, PizzaChili_Text());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> > const, PizzaChili_Text>::Type
getFibre(Index<TText, PizzaChili<TSpec> > const& me, PizzaChili_Text const) {
SEQAN_CHECKPOINT
    return value(me.text);
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> >, PizzaChili_Text>::Type
getFibre(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Text const) {
SEQAN_CHECKPOINT
    return value(me.text);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline bool
indexSupplied(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Compressed const) {
SEQAN_CHECKPOINT
    return me.index_handle != 0;
}

template <typename TText, typename TSpec>
inline bool
indexSupplied(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Text const) {
SEQAN_CHECKPOINT
    return length(value(me.text)) > 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline bool
indexSolveDependencies(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Compressed const) {
SEQAN_CHECKPOINT
    return indexSupplied(me, PizzaChili_Text());
}

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TText, typename TSpec>
    inline bool
    createPizzaChiliIndex(
        Index<TText, PizzaChili<TSpec> >& me,
        uchar_t* textstart,
        ulong_t textlength
    ) {
        typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
        // Read-only access, therefore safe cast.
        char* options = const_cast<char*>(impl::getOptionsString(me));
        impl::error_t e =
            TCodeProvider::build_index(textstart, textlength, options, &me.index_handle);

        if (e != 0) {
            SEQAN_REPORT(TCodeProvider::error_index(e));
            SEQAN_REPORT(options);
            me.index_handle = 0;
            return false;
        }

        value(me.text) = String<typename Value<TText>::Type, PizzaChili<TSpec> >(me.index_handle);

        return true;
    }
} // namespace impl

template <typename TText, typename TSpec>
inline bool
indexCreate(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Compressed const) {
SEQAN_CHECKPOINT
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
    typedef
        typename _RemoveConst<
            typename Index<TText, PizzaChili<TSpec> >::TValue
        >::Type alph_t;

    SEQAN_ASSERT(sizeof (alph_t) == 1);
    SEQAN_ASSERT((TYPECMP<typename IsSimple<alph_t>::Type, True>::VALUE));

    impl::clearIndex(me);

    impl::uchar_t* textstart =
        reinterpret_cast<impl::uchar_t*>(
            const_cast<alph_t*>(indexText(me).data_begin)
        );
    impl::ulong_t textlength = length(indexText(me));
    return impl::createPizzaChiliIndex(me, textstart, textlength);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec, typename TOtherText>
inline void
setIndexText(Index<TText, PizzaChili<TSpec> >& me, TOtherText& text) {
SEQAN_CHECKPOINT
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
    typedef
        typename _RemoveConst<
            typename Value<TOtherText>::Type
        >::Type alph_t;
    clear(me);

    /*
    SEQAN_ASSERT(IsContiguous<TOtherText>::VALUE)
    SEQAN_ASSERT(BitsPerValue<alph_t>::VALUE == 8)
    SEQAN_ASSERT((TYPECMP<typename IsSimple<alph_t>::Type, True>::VALUE));

    String<alph_t, CStyle> cstr = text;
    impl::uchar_t* textstart =
        reinterpret_cast<impl::uchar_t*>(
            const_cast<alph_t*>(static_cast<alph_t const*>(cstr))
        );
    impl::ulong_t textlength = length(text);
    impl::createPizzaChiliIndex(me, textstart, textlength);
    */
    getFibre(me, PizzaChili_Text()) = text;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline bool open(
    Index<TText, PizzaChili<TSpec> >& me,
    char const* filename
) {
SEQAN_CHECKPOINT
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
    clear(me);
    impl::error_t e =
        TCodeProvider::load_index(const_cast<char*>(filename), &me.index_handle);
    if (e != 0) {
        SEQAN_REPORT(TCodeProvider::error_index(e));
    }
    else
        value(me.text) = String<typename Value<TText>::Type, PizzaChili<TSpec> >(me.index_handle);
    return e == 0;
}

template <typename TText, typename TSpec>
inline bool save(
    Index<TText, PizzaChili<TSpec> >& me,
    char const* filename
) {
SEQAN_CHECKPOINT
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;
    // Before we can save the index we have to construct it, in case it isn't
    // already constructed.
    indexRequire(me, PizzaChili_Compressed());
    impl::error_t e =
        TCodeProvider::save_index(me.index_handle, const_cast<char*>(filename));
    if (e != 0) {
        SEQAN_REPORT(TCodeProvider::error_index(e));
    }
    return e == 0;
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PIZZACHILI_H
