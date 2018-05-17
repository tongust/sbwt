# Todo lists
1. [Build Index](#build-index)
2. [Support for Multiple Choromosomes Searching](#support-for-multiple)
3. [Hasing to Power The Seeding Phase](#hasing-to-power-the-seeding)
4. [Sort Index in Reverse-Complementary Way](#sort-index-in-reverse-complemetary-way)
5. [Load Index](#load-index)
6. [Reduce Memory Footprint in FM-index](#reduce-memory-footprint-in-fm-index)

## Build Index from Files
Should build second index from index files those built already instead of running on the fly.

##  Support for Multiple Choromosomes Searching
Currently, SBWT could not handle the multiple choromosomes sequence such human space geomoe.
To deal with that, **range search algorithm** should be implemented.

## Hashing to Power The Seeding Phase
Hashing up to 10 base pair to initialize the seeding.

## Sort Index in Reverse-Complementary Way
The second index would not used in reverse-complementary verification.

## Load index
Use memory mapping to load index instead of read word-wise.

## Reduce Memory Footprint in FM-index