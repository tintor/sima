#pragma once
#include <core/auto.h>
#include <core/murmur3.h>
#include <algorithm>
#include <cassert>
#include <string>
#include <string_view>
#include <cstring>

// clang-format off
struct prime_number_hash_policy {
		static size_t mod0(size_t) { return 0llu; }
		static size_t mod2(size_t hash) { return hash % 2llu; }
		static size_t mod3(size_t hash) { return hash % 3llu; }
		static size_t mod5(size_t hash) { return hash % 5llu; }
		static size_t mod7(size_t hash) { return hash % 7llu; }
		static size_t mod11(size_t hash) { return hash % 11llu; }
		static size_t mod13(size_t hash) { return hash % 13llu; }
		static size_t mod17(size_t hash) { return hash % 17llu; }
		static size_t mod23(size_t hash) { return hash % 23llu; }
		static size_t mod29(size_t hash) { return hash % 29llu; }
		static size_t mod37(size_t hash) { return hash % 37llu; }
		static size_t mod47(size_t hash) { return hash % 47llu; }
		static size_t mod59(size_t hash) { return hash % 59llu; }
		static size_t mod73(size_t hash) { return hash % 73llu; }
		static size_t mod97(size_t hash) { return hash % 97llu; }
		static size_t mod127(size_t hash) { return hash % 127llu; }
		static size_t mod151(size_t hash) { return hash % 151llu; }
		static size_t mod197(size_t hash) { return hash % 197llu; }
		static size_t mod251(size_t hash) { return hash % 251llu; }
		static size_t mod313(size_t hash) { return hash % 313llu; }
		static size_t mod397(size_t hash) { return hash % 397llu; }
		static size_t mod499(size_t hash) { return hash % 499llu; }
		static size_t mod631(size_t hash) { return hash % 631llu; }
		static size_t mod797(size_t hash) { return hash % 797llu; }
		static size_t mod1009(size_t hash) { return hash % 1009llu; }
		static size_t mod1259(size_t hash) { return hash % 1259llu; }
		static size_t mod1597(size_t hash) { return hash % 1597llu; }
		static size_t mod2011(size_t hash) { return hash % 2011llu; }
		static size_t mod2539(size_t hash) { return hash % 2539llu; }
		static size_t mod3203(size_t hash) { return hash % 3203llu; }
		static size_t mod4027(size_t hash) { return hash % 4027llu; }
		static size_t mod5087(size_t hash) { return hash % 5087llu; }
		static size_t mod6421(size_t hash) { return hash % 6421llu; }
		static size_t mod8089(size_t hash) { return hash % 8089llu; }
		static size_t mod10193(size_t hash) { return hash % 10193llu; }
		static size_t mod12853(size_t hash) { return hash % 12853llu; }
		static size_t mod16193(size_t hash) { return hash % 16193llu; }
		static size_t mod20399(size_t hash) { return hash % 20399llu; }
		static size_t mod25717(size_t hash) { return hash % 25717llu; }
		static size_t mod32401(size_t hash) { return hash % 32401llu; }
		static size_t mod40823(size_t hash) { return hash % 40823llu; }
		static size_t mod51437(size_t hash) { return hash % 51437llu; }
		static size_t mod64811(size_t hash) { return hash % 64811llu; }
		static size_t mod81649(size_t hash) { return hash % 81649llu; }
		static size_t mod102877(size_t hash) { return hash % 102877llu; }
		static size_t mod129607(size_t hash) { return hash % 129607llu; }
		static size_t mod163307(size_t hash) { return hash % 163307llu; }
		static size_t mod205759(size_t hash) { return hash % 205759llu; }
		static size_t mod259229(size_t hash) { return hash % 259229llu; }
		static size_t mod326617(size_t hash) { return hash % 326617llu; }
		static size_t mod411527(size_t hash) { return hash % 411527llu; }
		static size_t mod518509(size_t hash) { return hash % 518509llu; }
		static size_t mod653267(size_t hash) { return hash % 653267llu; }
		static size_t mod823117(size_t hash) { return hash % 823117llu; }
		static size_t mod1037059(size_t hash) { return hash % 1037059llu; }
		static size_t mod1306601(size_t hash) { return hash % 1306601llu; }
		static size_t mod1646237(size_t hash) { return hash % 1646237llu; }
		static size_t mod2074129(size_t hash) { return hash % 2074129llu; }
		static size_t mod2613229(size_t hash) { return hash % 2613229llu; }
		static size_t mod3292489(size_t hash) { return hash % 3292489llu; }
		static size_t mod4148279(size_t hash) { return hash % 4148279llu; }
		static size_t mod5226491(size_t hash) { return hash % 5226491llu; }
		static size_t mod6584983(size_t hash) { return hash % 6584983llu; }
		static size_t mod8296553(size_t hash) { return hash % 8296553llu; }
		static size_t mod10453007(size_t hash) { return hash % 10453007llu; }
		static size_t mod13169977(size_t hash) { return hash % 13169977llu; }
		static size_t mod16593127(size_t hash) { return hash % 16593127llu; }
		static size_t mod20906033(size_t hash) { return hash % 20906033llu; }
		static size_t mod26339969(size_t hash) { return hash % 26339969llu; }
		static size_t mod33186281(size_t hash) { return hash % 33186281llu; }
		static size_t mod41812097(size_t hash) { return hash % 41812097llu; }
		static size_t mod52679969(size_t hash) { return hash % 52679969llu; }
		static size_t mod66372617(size_t hash) { return hash % 66372617llu; }
		static size_t mod83624237(size_t hash) { return hash % 83624237llu; }
		static size_t mod105359939(size_t hash) { return hash % 105359939llu; }
		static size_t mod132745199(size_t hash) { return hash % 132745199llu; }
		static size_t mod167248483(size_t hash) { return hash % 167248483llu; }
		static size_t mod210719881(size_t hash) { return hash % 210719881llu; }
		static size_t mod265490441(size_t hash) { return hash % 265490441llu; }
		static size_t mod334496971(size_t hash) { return hash % 334496971llu; }
		static size_t mod421439783(size_t hash) { return hash % 421439783llu; }
		static size_t mod530980861(size_t hash) { return hash % 530980861llu; }
		static size_t mod668993977(size_t hash) { return hash % 668993977llu; }
		static size_t mod842879579(size_t hash) { return hash % 842879579llu; }
		static size_t mod1061961721(size_t hash) { return hash % 1061961721llu; }
		static size_t mod1337987929(size_t hash) { return hash % 1337987929llu; }
		static size_t mod1685759167(size_t hash) { return hash % 1685759167llu; }
		static size_t mod2123923447(size_t hash) { return hash % 2123923447llu; }
		static size_t mod2675975881(size_t hash) { return hash % 2675975881llu; }
		static size_t mod3371518343(size_t hash) { return hash % 3371518343llu; }
		static size_t mod4247846927(size_t hash) { return hash % 4247846927llu; }
		static size_t mod5351951779(size_t hash) { return hash % 5351951779llu; }
		static size_t mod6743036717(size_t hash) { return hash % 6743036717llu; }
		static size_t mod8495693897(size_t hash) { return hash % 8495693897llu; }
		static size_t mod10703903591(size_t hash) { return hash % 10703903591llu; }
		static size_t mod13486073473(size_t hash) { return hash % 13486073473llu; }
		static size_t mod16991387857(size_t hash) { return hash % 16991387857llu; }
		static size_t mod21407807219(size_t hash) { return hash % 21407807219llu; }
		static size_t mod26972146961(size_t hash) { return hash % 26972146961llu; }
		static size_t mod33982775741(size_t hash) { return hash % 33982775741llu; }
		static size_t mod42815614441(size_t hash) { return hash % 42815614441llu; }
		static size_t mod53944293929(size_t hash) { return hash % 53944293929llu; }
		static size_t mod67965551447(size_t hash) { return hash % 67965551447llu; }
		static size_t mod85631228929(size_t hash) { return hash % 85631228929llu; }
		static size_t mod107888587883(size_t hash) { return hash % 107888587883llu; }
		static size_t mod135931102921(size_t hash) { return hash % 135931102921llu; }
		static size_t mod171262457903(size_t hash) { return hash % 171262457903llu; }
		static size_t mod215777175787(size_t hash) { return hash % 215777175787llu; }
		static size_t mod271862205833(size_t hash) { return hash % 271862205833llu; }
		static size_t mod342524915839(size_t hash) { return hash % 342524915839llu; }
		static size_t mod431554351609(size_t hash) { return hash % 431554351609llu; }
		static size_t mod543724411781(size_t hash) { return hash % 543724411781llu; }
		static size_t mod685049831731(size_t hash) { return hash % 685049831731llu; }
		static size_t mod863108703229(size_t hash) { return hash % 863108703229llu; }
		static size_t mod1087448823553(size_t hash) { return hash % 1087448823553llu; }
		static size_t mod1370099663459(size_t hash) { return hash % 1370099663459llu; }
		static size_t mod1726217406467(size_t hash) { return hash % 1726217406467llu; }
		static size_t mod2174897647073(size_t hash) { return hash % 2174897647073llu; }
		static size_t mod2740199326961(size_t hash) { return hash % 2740199326961llu; }
		static size_t mod3452434812973(size_t hash) { return hash % 3452434812973llu; }
		static size_t mod4349795294267(size_t hash) { return hash % 4349795294267llu; }
		static size_t mod5480398654009(size_t hash) { return hash % 5480398654009llu; }
		static size_t mod6904869625999(size_t hash) { return hash % 6904869625999llu; }
		static size_t mod8699590588571(size_t hash) { return hash % 8699590588571llu; }
		static size_t mod10960797308051(size_t hash) { return hash % 10960797308051llu; }
		static size_t mod13809739252051(size_t hash) { return hash % 13809739252051llu; }
		static size_t mod17399181177241(size_t hash) { return hash % 17399181177241llu; }
		static size_t mod21921594616111(size_t hash) { return hash % 21921594616111llu; }
		static size_t mod27619478504183(size_t hash) { return hash % 27619478504183llu; }
		static size_t mod34798362354533(size_t hash) { return hash % 34798362354533llu; }
		static size_t mod43843189232363(size_t hash) { return hash % 43843189232363llu; }
		static size_t mod55238957008387(size_t hash) { return hash % 55238957008387llu; }
		static size_t mod69596724709081(size_t hash) { return hash % 69596724709081llu; }
		static size_t mod87686378464759(size_t hash) { return hash % 87686378464759llu; }
		static size_t mod110477914016779(size_t hash) { return hash % 110477914016779llu; }
		static size_t mod139193449418173(size_t hash) { return hash % 139193449418173llu; }
		static size_t mod175372756929481(size_t hash) { return hash % 175372756929481llu; }
		static size_t mod220955828033581(size_t hash) { return hash % 220955828033581llu; }
		static size_t mod278386898836457(size_t hash) { return hash % 278386898836457llu; }
		static size_t mod350745513859007(size_t hash) { return hash % 350745513859007llu; }
		static size_t mod441911656067171(size_t hash) { return hash % 441911656067171llu; }
		static size_t mod556773797672909(size_t hash) { return hash % 556773797672909llu; }
		static size_t mod701491027718027(size_t hash) { return hash % 701491027718027llu; }
		static size_t mod883823312134381(size_t hash) { return hash % 883823312134381llu; }
		static size_t mod1113547595345903(size_t hash) { return hash % 1113547595345903llu; }
		static size_t mod1402982055436147(size_t hash) { return hash % 1402982055436147llu; }
		static size_t mod1767646624268779(size_t hash) { return hash % 1767646624268779llu; }
		static size_t mod2227095190691797(size_t hash) { return hash % 2227095190691797llu; }
		static size_t mod2805964110872297(size_t hash) { return hash % 2805964110872297llu; }
		static size_t mod3535293248537579(size_t hash) { return hash % 3535293248537579llu; }
		static size_t mod4454190381383713(size_t hash) { return hash % 4454190381383713llu; }
		static size_t mod5611928221744609(size_t hash) { return hash % 5611928221744609llu; }
		static size_t mod7070586497075177(size_t hash) { return hash % 7070586497075177llu; }
		static size_t mod8908380762767489(size_t hash) { return hash % 8908380762767489llu; }
		static size_t mod11223856443489329(size_t hash) { return hash % 11223856443489329llu; }
		static size_t mod14141172994150357(size_t hash) { return hash % 14141172994150357llu; }
		static size_t mod17816761525534927(size_t hash) { return hash % 17816761525534927llu; }
		static size_t mod22447712886978529(size_t hash) { return hash % 22447712886978529llu; }
		static size_t mod28282345988300791(size_t hash) { return hash % 28282345988300791llu; }
		static size_t mod35633523051069991(size_t hash) { return hash % 35633523051069991llu; }
		static size_t mod44895425773957261(size_t hash) { return hash % 44895425773957261llu; }
		static size_t mod56564691976601587(size_t hash) { return hash % 56564691976601587llu; }
		static size_t mod71267046102139967(size_t hash) { return hash % 71267046102139967llu; }
		static size_t mod89790851547914507(size_t hash) { return hash % 89790851547914507llu; }
		static size_t mod113129383953203213(size_t hash) { return hash % 113129383953203213llu; }
		static size_t mod142534092204280003(size_t hash) { return hash % 142534092204280003llu; }
		static size_t mod179581703095829107(size_t hash) { return hash % 179581703095829107llu; }
		static size_t mod226258767906406483(size_t hash) { return hash % 226258767906406483llu; }
		static size_t mod285068184408560057(size_t hash) { return hash % 285068184408560057llu; }
		static size_t mod359163406191658253(size_t hash) { return hash % 359163406191658253llu; }
		static size_t mod452517535812813007(size_t hash) { return hash % 452517535812813007llu; }
		static size_t mod570136368817120201(size_t hash) { return hash % 570136368817120201llu; }
		static size_t mod718326812383316683(size_t hash) { return hash % 718326812383316683llu; }
		static size_t mod905035071625626043(size_t hash) { return hash % 905035071625626043llu; }
		static size_t mod1140272737634240411(size_t hash) { return hash % 1140272737634240411llu; }
		static size_t mod1436653624766633509(size_t hash) { return hash % 1436653624766633509llu; }
		static size_t mod1810070143251252131(size_t hash) { return hash % 1810070143251252131llu; }
		static size_t mod2280545475268481167(size_t hash) { return hash % 2280545475268481167llu; }
		static size_t mod2873307249533267101(size_t hash) { return hash % 2873307249533267101llu; }
		static size_t mod3620140286502504283(size_t hash) { return hash % 3620140286502504283llu; }
		static size_t mod4561090950536962147(size_t hash) { return hash % 4561090950536962147llu; }
		static size_t mod5746614499066534157(size_t hash) { return hash % 5746614499066534157llu; }
		static size_t mod7240280573005008577(size_t hash) { return hash % 7240280573005008577llu; }
		static size_t mod9122181901073924329(size_t hash) { return hash % 9122181901073924329llu; }
		static size_t mod11493228998133068689(size_t hash) { return hash % 11493228998133068689llu; }
		static size_t mod14480561146010017169(size_t hash) { return hash % 14480561146010017169llu; }
		static size_t mod18446744073709551557(size_t hash) { return hash % 18446744073709551557llu; }

		using mod_function = size_t (*)(size_t);

		static constexpr const size_t prime_list[] = {
				2llu, 3llu, 5llu, 7llu, 11llu, 13llu, 17llu, 23llu, 29llu, 37llu, 47llu,
				59llu, 73llu, 97llu, 127llu, 151llu, 197llu, 251llu, 313llu, 397llu,
				499llu, 631llu, 797llu, 1009llu, 1259llu, 1597llu, 2011llu, 2539llu,
				3203llu, 4027llu, 5087llu, 6421llu, 8089llu, 10193llu, 12853llu, 16193llu,
				20399llu, 25717llu, 32401llu, 40823llu, 51437llu, 64811llu, 81649llu,
				102877llu, 129607llu, 163307llu, 205759llu, 259229llu, 326617llu,
				411527llu, 518509llu, 653267llu, 823117llu, 1037059llu, 1306601llu,
				1646237llu, 2074129llu, 2613229llu, 3292489llu, 4148279llu, 5226491llu,
				6584983llu, 8296553llu, 10453007llu, 13169977llu, 16593127llu, 20906033llu,
				26339969llu, 33186281llu, 41812097llu, 52679969llu, 66372617llu,
				83624237llu, 105359939llu, 132745199llu, 167248483llu, 210719881llu,
				265490441llu, 334496971llu, 421439783llu, 530980861llu, 668993977llu,
				842879579llu, 1061961721llu, 1337987929llu, 1685759167llu, 2123923447llu,
				2675975881llu, 3371518343llu, 4247846927llu, 5351951779llu, 6743036717llu,
				8495693897llu, 10703903591llu, 13486073473llu, 16991387857llu,
				21407807219llu, 26972146961llu, 33982775741llu, 42815614441llu,
				53944293929llu, 67965551447llu, 85631228929llu, 107888587883llu,
				135931102921llu, 171262457903llu, 215777175787llu, 271862205833llu,
				342524915839llu, 431554351609llu, 543724411781llu, 685049831731llu,
				863108703229llu, 1087448823553llu, 1370099663459llu, 1726217406467llu,
				2174897647073llu, 2740199326961llu, 3452434812973llu, 4349795294267llu,
				5480398654009llu, 6904869625999llu, 8699590588571llu, 10960797308051llu,
				13809739252051llu, 17399181177241llu, 21921594616111llu, 27619478504183llu,
				34798362354533llu, 43843189232363llu, 55238957008387llu, 69596724709081llu,
				87686378464759llu, 110477914016779llu, 139193449418173llu,
				175372756929481llu, 220955828033581llu, 278386898836457llu,
				350745513859007llu, 441911656067171llu, 556773797672909llu,
				701491027718027llu, 883823312134381llu, 1113547595345903llu,
				1402982055436147llu, 1767646624268779llu, 2227095190691797llu,
				2805964110872297llu, 3535293248537579llu, 4454190381383713llu,
				5611928221744609llu, 7070586497075177llu, 8908380762767489llu,
				11223856443489329llu, 14141172994150357llu, 17816761525534927llu,
				22447712886978529llu, 28282345988300791llu, 35633523051069991llu,
				44895425773957261llu, 56564691976601587llu, 71267046102139967llu,
				89790851547914507llu, 113129383953203213llu, 142534092204280003llu,
				179581703095829107llu, 226258767906406483llu, 285068184408560057llu,
				359163406191658253llu, 452517535812813007llu, 570136368817120201llu,
				718326812383316683llu, 905035071625626043llu, 1140272737634240411llu,
				1436653624766633509llu, 1810070143251252131llu, 2280545475268481167llu,
				2873307249533267101llu, 3620140286502504283llu, 4561090950536962147llu,
				5746614499066534157llu, 7240280573005008577llu, 9122181901073924329llu,
				11493228998133068689llu, 14480561146010017169llu, 18446744073709551557llu
		};

		mod_function next_size_over(size_t & size) const {
				// prime numbers generated by the following method:
				// 1. start with a prime p = 2
				// 2. go to wolfram alpha and get p = NextPrime(2 * p)
				// 3. repeat 2. until you overflow 64 bits
				// you now have large gaps which you would hit if somebody called reserve() with an unlucky number.
				// 4. to fill the gaps for every prime p go to wolfram alpha and get ClosestPrime(p * 2^(1/3)) and ClosestPrime(p * 2^(2/3)) and put those in the gaps
				// 5. get PrevPrime(2^64) and put it at the end
				static constexpr size_t (* const mod_functions[])(size_t) = {
						&mod0, &mod2, &mod3, &mod5, &mod7, &mod11, &mod13, &mod17, &mod23, &mod29, &mod37,
						&mod47, &mod59, &mod73, &mod97, &mod127, &mod151, &mod197, &mod251, &mod313, &mod397,
						&mod499, &mod631, &mod797, &mod1009, &mod1259, &mod1597, &mod2011, &mod2539, &mod3203,
						&mod4027, &mod5087, &mod6421, &mod8089, &mod10193, &mod12853, &mod16193, &mod20399,
						&mod25717, &mod32401, &mod40823, &mod51437, &mod64811, &mod81649, &mod102877,
						&mod129607, &mod163307, &mod205759, &mod259229, &mod326617, &mod411527, &mod518509,
						&mod653267, &mod823117, &mod1037059, &mod1306601, &mod1646237, &mod2074129,
						&mod2613229, &mod3292489, &mod4148279, &mod5226491, &mod6584983, &mod8296553,
						&mod10453007, &mod13169977, &mod16593127, &mod20906033, &mod26339969, &mod33186281,
						&mod41812097, &mod52679969, &mod66372617, &mod83624237, &mod105359939, &mod132745199,
						&mod167248483, &mod210719881, &mod265490441, &mod334496971, &mod421439783,
						&mod530980861, &mod668993977, &mod842879579, &mod1061961721, &mod1337987929,
						&mod1685759167, &mod2123923447, &mod2675975881, &mod3371518343, &mod4247846927,
						&mod5351951779, &mod6743036717, &mod8495693897, &mod10703903591, &mod13486073473,
						&mod16991387857, &mod21407807219, &mod26972146961, &mod33982775741, &mod42815614441,
						&mod53944293929, &mod67965551447, &mod85631228929, &mod107888587883, &mod135931102921,
						&mod171262457903, &mod215777175787, &mod271862205833, &mod342524915839,
						&mod431554351609, &mod543724411781, &mod685049831731, &mod863108703229,
						&mod1087448823553, &mod1370099663459, &mod1726217406467, &mod2174897647073,
						&mod2740199326961, &mod3452434812973, &mod4349795294267, &mod5480398654009,
						&mod6904869625999, &mod8699590588571, &mod10960797308051, &mod13809739252051,
						&mod17399181177241, &mod21921594616111, &mod27619478504183, &mod34798362354533,
						&mod43843189232363, &mod55238957008387, &mod69596724709081, &mod87686378464759,
						&mod110477914016779, &mod139193449418173, &mod175372756929481, &mod220955828033581,
						&mod278386898836457, &mod350745513859007, &mod441911656067171, &mod556773797672909,
						&mod701491027718027, &mod883823312134381, &mod1113547595345903, &mod1402982055436147,
						&mod1767646624268779, &mod2227095190691797, &mod2805964110872297, &mod3535293248537579,
						&mod4454190381383713, &mod5611928221744609, &mod7070586497075177, &mod8908380762767489,
						&mod11223856443489329, &mod14141172994150357, &mod17816761525534927,
						&mod22447712886978529, &mod28282345988300791, &mod35633523051069991,
						&mod44895425773957261, &mod56564691976601587, &mod71267046102139967,
						&mod89790851547914507, &mod113129383953203213, &mod142534092204280003,
						&mod179581703095829107, &mod226258767906406483, &mod285068184408560057,
						&mod359163406191658253, &mod452517535812813007, &mod570136368817120201,
						&mod718326812383316683, &mod905035071625626043, &mod1140272737634240411,
						&mod1436653624766633509, &mod1810070143251252131, &mod2280545475268481167,
						&mod2873307249533267101, &mod3620140286502504283, &mod4561090950536962147,
						&mod5746614499066534157, &mod7240280573005008577, &mod9122181901073924329,
						&mod11493228998133068689, &mod14480561146010017169, &mod18446744073709551557
				};
				const size_t * found = std::lower_bound(std::begin(prime_list), std::end(prime_list) - 1, size);
				size = *found;
				return mod_functions[1 + found - prime_list];
		}

		size_t capacity() const {
				return 0; // TODO
		}

		mod_function fn = &mod0;
};
// clang-format on

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

inline bool aligned(const void* p, size_t a) {
	return reinterpret_cast<size_t>(p) % a == 0;
}

#ifdef __clang__
constexpr bool _clang_ = true;
#else
constexpr bool _clang_ = false;
#endif

struct fast_string_view {
	constexpr static bool special = true;

	fast_string_view(const char* s) {
		uint32_t l = strlen(s);
		init(l, s, l + 1);
	}

	fast_string_view(const std::string& s) {
		// clang std::string stores \0 char after the last reserved byte
		init(s.size(), s.data(), s.capacity() + (_clang_ ? 1 : 0));
	}

	fast_string_view(std::string_view s) { init(s.size(), s.data(), s.size()); }

	void init(uint32_t len, const char* data, uint32_t capacity) {
		assert(data != nullptr);
		ON_SCOPE_EXIT(assert(hash(0) == MurmurHash3_x64_128(data, len, 0)));
		_len = len;
		_data = reinterpret_cast<const uint64_t*>(data);
		if (_len % 8 == 0) {
			if (_len == 0) {
				_b = 0;
				_a = 0;
				return;
			}
			if (_len == 8) {
				_b = 0;
				_a = _data[0];
				return;
			}
			if (special && _len == 16) {
				_b = _data[0];
				_a = _data[1];
				_data = nullptr;
				return;
			}
			_b = _len / 8 - 1;
			_a = _data[_b];
			return;
		}
		if (aligned(_data, 8)) {
			_b = _len / 8;
			if (_b * 8 + 8 <= capacity) {
				_a = _data[_b]; // safe to read past the end of
				// string
				if (special && _len < 16) {
					_b = _data[0];
					_data = nullptr;
				}
				return;
			}
			// fallthrough
		}
		if (_len < 8) {
			_b = 0;
			_a = load64(data, _len);
			return;
		}
		if (special && _len < 16) {
			_b = *_data;
			_a = load64(data + 8, _len);
			_data = nullptr;
			return;
		}
		_b = _len / 8;
		_a = load64(data + _b * 8, _len % 8);
	}

	// safely load first len bytes from s into uint64_t (pad remaining bytes
	// with zeroes)
	static uint64_t load64(const char* s, int len) {
		assert(1 <= len && len <= 7);
		uint64_t k = 0;
		switch (len) {
		case 7:
			k ^= ((uint64_t)s[6]) << 48;
		case 6:
			k ^= ((uint64_t)s[5]) << 40;
		case 5:
			k ^= ((uint64_t)s[4]) << 32;
		case 4:
			return k ^ *reinterpret_cast<const uint32_t*>(s);
		case 3:
			k ^= ((uint64_t)s[2]) << 16;
		case 2:
			k ^= ((uint64_t)s[1]) << 8;
		case 1:
			k ^= ((uint64_t)s[0]) << 0;
		};
		return k;
	}

	fast_string_view& set_hash(size_t hash) {
		_hash = hash;
		return *this;
	}

	size_t hash(uint32_t seed) {
		if (_hash != 0) {
			return _hash;
		}
		Murmur3_x64_128 murmur;
		murmur.reset(seed);
		if (special && _data == nullptr) {
			murmur.tail(_b, _a);
		} else {
			for (uint32_t i = 0; i < _b / 2 * 2; i += 2)
				murmur.append(_data[i], _data[i + 1]);
			if (_b & 1)
				murmur.tail(_data[_b - 1], _a);
			else
				murmur.tail(_a, 0);
		}
		murmur.finalize(_len);
		return murmur.h1 ^ murmur.h2;
	}

	// key is 8 byte aligned and is allocated with size which is multiple of
	// 8 bytes
	bool equals(const uint64_t* key, uint32_t key_len) {
		if (_len != key_len)
			return false;
		if (special && _data == nullptr)
			return _b == key[0] && _a == key[1];
		switch (_b) {
		case 0:
			return _a == key[0];
		case 1:
			return _data[0] == key[0] && _a == key[1];
		case 2:
			return _data[0] == key[0] && _data[1] == key[1] && _a == key[2];
		case 3:
			return _data[0] == key[0] && _data[1] == key[1] &&
				   _data[2] == key[2] && _a == key[3];
		case 4:
			return _data[0] == key[0] && _data[1] == key[1] &&
				   _data[2] == key[2] && _data[3] == key[3] && _a == key[4];
		}
		for (uint32_t i = 0; i < _b; i++)
			if (_data[i] != key[i])
				return false;
		return _a == key[_b];
	}

	uint32_t size() { return _len; }

	void write(uint64_t* key) {
		if (special && _data == nullptr) {
			key[0] = _b;
			key[1] = _a;
			return;
		}
		for (uint32_t i = 0; i < _b; i++)
			key[0] = _data[i];
		key[_b] = _a;
	}

  private:
	const uint64_t* _data;
	uint64_t _a, _b;
	// without special replace _b with uint32_t _blocks
	uint32_t _len;
	size_t _hash = 0;
};

// similar to std::unordered_map<std::string, V>, but uses less memory
// string keys are inlined into nodes
// uses murmur3_x64_128 hash
//
// TODO implement std::swap
// TODO printing of container to ostream
// TODO when removing elements, save node allocations for future inserts (of the
// same key size)
template <typename V, int milli_load_factor = 1000, uint32_t hash_seed = 0>
class string_node_hash_map {
	constexpr static float load_factor = milli_load_factor * 1e-3;

  public:
	struct node {
		node* _next;
		V _value;
		uint32_t _key_size;
		uint64_t _key[1] __attribute__((aligned(8)));

		node(node* next, fast_string_view s, const V& value)
			: _next(next), _value(value) {
			_key_size = s.size();
			s.write(_key);
		}
		static void* operator new(size_t sz, size_t key_size) {
			key_size = (std::max<size_t>(1, key_size) + 7) / 8 * 8;
			return ::operator new(sz + key_size);
		}
		std::string_view key() const {
			return std::string_view(reinterpret_cast<const char*>(_key),
									_key_size);
		}
		size_t hash() {
			/*Murmur3_x64_128 murmur;
			murmur.reset(hash_seed);
			for (uint32_t i = 0; i < _key_size / 8; i += 2)
				murmur.append(_key[i], _key[i + 1]);
			if (_key_size % 16 <= 8)
				murmur.tail(_key[last], _a);
			murmur.finalize(_key_size);
			return murmur.h1 ^ murmur.h2;*/
			// TODO
			return 0;
		}
	};

	class iterator {
	  public:
		iterator() {}
		iterator(node** p, node** e, node* n) : _p(p), _e(e), _n(n) {}

		std::pair<std::string_view, V&> operator*() {
			return {_n->key(), _n->_value};
		}
		std::pair<std::string_view, const V&> operator*() const {
			return {_n->key(), _n->_value};
		}

		std::string_view key() const { return _n->key(); }
		const V& value() const { return _n->_value; }
		V& value() { return _n->_value; }

		void operator++() {
			_n = _n->_next;
			if (_n == nullptr) {
				_p += 1;
				if (_p == _e)
					return;
				_n = *_p;
			}
		}
		void operator++(int) { operator++(); }
		bool operator!=(const iterator& o) {
			return _p != o._p || _e != o._e || _n != o._n;
		}

	  private:
		node** _p;
		node** _e;
		node* _n;
	};

	class const_iterator {
	  public:
		const_iterator() {}
		const_iterator(const node** p, const node** e, const node* n)
			: _p(p), _e(e), _n(n) {}

		std::pair<std::string_view, const V&> operator*() const {
			return {_n->key(), _n->_value};
		}

		std::string_view key() const { return _n->key(); }
		const V& value() const { return _n->_value; }

		void operator++() {
			_n = _n->next;
			if (_n == nullptr) {
				_p += 1;
				if (_p == _e)
					return;
				_n = *_p;
			}
		}
		void operator++(int) { operator++(); }
		bool operator!=(const iterator& o) {
			return _p != o._p || _e != o._e || _n != o._n;
		}

	  private:
		const node** _p;
		const node** _e;
		const node* _n;
	};

	string_node_hash_map() : _size(0), _entries(nullptr) {}

	// TODO copy constructor
	// TODO move constructor
	// TODO assignment operator
	// TODO move assignment operator

	~string_node_hash_map() {
		if (_entries) {
			size_t capacity = _hash_policy.capacity();
			for (size_t i = 0; i < capacity; i++) {
				node* e = _entries[i];
				while (e) {
					node* n = e->_next;
					try {
						delete e;
					} catch (...) {
						// swallow exception to avoid
						// leaking memory
					}
					e = n;
				}
			}
			std::free(_entries);
		}
	}

	const_iterator begin() const {
		return const_iterator(_entries, _entries + capacity(), nullptr);
	}
	const_iterator end() const {
		return const_iterator(_entries + capacity(), _entries + capacity(),
							  nullptr);
	}
	iterator begin() {
		return iterator(_entries, _entries + capacity(), nullptr);
	}
	iterator end() {
		return iterator(_entries + capacity(), _entries + capacity(), nullptr);
	}

	size_t size() const { return size; }
	size_t capacity() const { return _hash_policy.capacity(); }

	void reserve(size_t new_capacity) {
		size_t capacity = _hash_policy.capacity();
		if (capacity == 0) {
			shrink_to_fit();
			return;
		}
		auto new_mod_fn = _hash_policy.next_size_over(/*ref*/ new_capacity);
		resize(capacity, new_capacity, new_mod_fn);
	}

	void clear() {
		size_t capacity = _hash_policy.capacity();
		for (size_t i = 0; i < capacity; i++) {
			while (_entries[i]) {
				node* n = _entries[i]->_next;
				delete _entries[i];
				_size -= 1;
				_entries[i] = n;
			}
		}
		assert(_size == 0);
	}

	void shrink_to_fit() {
		if (_size == 0) {
			_hash_policy.fn = &prime_number_hash_policy::mod0;
			std::free(_entries);
			_entries = nullptr;
			return;
		}

		size_t capacity = _hash_policy.capacity();
		size_t new_capacity = _size / load_factor;
		auto new_mod_fn = _hash_policy.next_size_over(/*ref*/ new_capacity);
		if (new_capacity != capacity) {
			resize(capacity, new_capacity, new_mod_fn);
		}
	}

	void insert(fast_string_view key, V value) {
		if (unlikely(_entries == nullptr)) {
			size_t new_capacity = 4;
			_hash_policy.fn = _hash_policy.next_size_over(/*ref*/ new_capacity);
			_entries = reinterpret_cast<node**>(
				std::calloc(new_capacity, sizeof(node*)));
			if (_entries == nullptr)
				throw std::bad_alloc();
		}

		size_t h = index(key.hash(hash_seed));
		node* e = _entries[h];
		for (node* n = e; n; n = n->_next)
			if (key.equals(n->_key, n->_key_size)) {
				n->_value = value;
				return;
			}
		if (_size + 1 >= load_factor * _hash_policy.capacity())
			grow();
		_entries[h] = new (key.size()) node(e, key, value);
		assert(_entries[h]->hash() == key.hash(hash_seed));
		// increment size after allocation in case of exception
		_size += 1;
	}

	bool remove(fast_string_view key) {
		size_t h = index(key.hash(hash_seed));
		node* e = _entries[h];
		if (e == nullptr)
			return false;
		if (key.equals(e->_key, e->_key_size)) {
			_entries[h] = e->next;
			delete e;
			_size -= 1;
			return true;
		}
		for (node* n = e->_next; n; e = n, n = n->_next)
			if (key.equals(n->_key, n->_key_size)) {
				e->_next = n->_next;
				delete n;
				_size -= 1;
				return true;
			}
		return false;
	}

	// bool contains(const std::string& key, size_t hash)  {
	bool contains(fast_string_view s) {
		size_t h = index(s.hash(hash_seed));
		for (node* n = _entries[h]; n; n = n->_next)
			if (s.equals(n->_key, n->_key_size))
				return true;
		return false;
	}

	const V& at(fast_string_view key) const {
		size_t h = index(key.hash(hash_seed));
		for (node* n = _entries[h]; n; n = n->_next)
			if (key.equals(n->_key, n->_key_size))
				return n->_value;
		throw std::out_of_range("key");
	}

	V& at(fast_string_view key) {
		size_t h = index(key.hash(hash_seed));
		for (node* n = _entries[h]; n; n = n->_next)
			if (key.equals(n->_key, n->_key_size))
				return n->_value;
		throw std::out_of_range("key");
	}

  private:
	void grow() {
		size_t capacity = _hash_policy.capacity();
		size_t new_capacity = capacity;
		auto new_mod_fn = _hash_policy.next_size_over(/*ref*/ new_capacity);
		resize(capacity, new_capacity, new_mod_fn);
	}

	void resize(size_t capacity, size_t new_capacity,
				prime_number_hash_policy::mod_function new_mod_fn) {
		node** new_entries =
			reinterpret_cast<node**>(std::calloc(new_capacity, sizeof(node*)));
		if (!new_entries)
			throw std::bad_alloc();
		// no exceptions after this line

		for (size_t i = 0; i < capacity; i++) {
			node* p = _entries[i];
			while (p) {
				node* n = p->_next;
				size_t h = new_mod_fn(p->hash());
				p->_next = new_entries[h];
				new_entries[h] = p;
				p = n;
			}
		}

		_hash_policy.fn = new_mod_fn;
		std::free(_entries);
		_entries = new_entries;
	}

	size_t index(size_t hash) { return _hash_policy.fn(hash); }

	size_t _size;
	node** _entries; // uses calloc/free
	prime_number_hash_policy _hash_policy;
};
