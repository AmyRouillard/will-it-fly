#ifndef __STYLE_MOM_HPP_INCLUDED__
#define __STYLE_MOM_HPP_INCLUDED__

// Will It Fly - coding guidelines, based on [SFT]tech coding guidelines
// ---------------------------------------------------------------------
//
// Please also have a look at `mom.cpp`!
// To make the code look nice and be readable:
//  * Sort things alphabetically.
//  * Align things beautifully!
//  * If you write documentation correctly, then:
//    * It's parsable by Doxygen!
//    * -> Cool generated docs!
//  * Look at the following example code.
//  * Look at the "mom.cpp" file for a implementation example.
//  * Comment guideline: Use punctuation!

#include <cstddef>
#include <cstdint>

// Try to keep the code < 120 chars in width (no hard limit...).
// All documentation is in the header.

// No swearing in comments!
// No preprocessor shit crap!
constexpr int DEFAULT_MOM_MASS = 9001;

// Doxygen will then generate nice shiny documentation files.
// Exported (API) documentation of all source code is in the header file:

/**
 * You may document a namespace if it's something special.
 * This documentation will also show up in generated docs then.
 */
namespace elts
{
namespace mom
{
// <- The namespace does not increase the indentation level!


/**
 * Class about your mother.
 * Manages all kinds of useful information about your mother.
 *
 * Classes are written like_this_c!
 */
class your_mom_c
{
public:
	/**
	 * Creates your mom.
	 * @param mass: Initial mass for this mother.
	 */
	your_mom_c(int mass = DEFAULT_MOM_MASS);

	~your_mom_c();
	// Note that the destructor has no doc,
	// Because that would be entirely unnecessary (in most cases).

	/**
	 * Prints your mom's current mass to stdout.
	 */
	void print_mass();

	/**
	 * Ticks your mom.
	 * Simulates her for one physics frame.
	 *
	 * @see mass
	 * @returns true if the simulation was successful.
	 */
	bool tick();

protected:
	// To have doxygen parse one-liners, write the following:
	int64_t mass; //!< Your mom's mass (solar masses).

	/**
	 * Your mom's number of non-artificial satellites.
	 * Scales well with your mom's mass.
	 * @see mass.
	 */
	int number_of_moons;

	/**
	 * Your mom's stupidity.
	 * Measured in unit of tolerating DRM.
	 */
	int stupidity;

	static int64_t sum_mom_masses; //!< Sum of all mom's masses in the world.
};


} // namespace elts
} // namespace mom


#endif // __STYLE_MOM_HPP_INCLUDED__
