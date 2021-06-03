#include <string.h>
#include "comb.h"
#include "curve.h"
#include "constant_time.h"
#include "f11_260.h"

// Comb set for base point. Used for fast signatures.
sabs_comb_set_t base_comb = {
  .combs = {
    {
      .table = {
        {
          .x = {
            0, 0x18c7a32, 0x335dfa6, 0x847151, 0x30e53d8, 0x3e91e82,
            0x8c8cf6, 0x1cf74cf, 0x1fa9f5e, 0x2ab354f, 0x19c9964, 0,
          },
          .dt = {
            0, 0x2266577, 0x58e652, 0xc88cc2, 0xd6883b, 0x1adbf4a,
            0x139a6da, 0x154d236, 0x26d98bc, 0x94fbf1, 0x31f79c1, 0,
          },
          .y = {
            0, 0x1ffe5ee, 0x257d19c, 0x2940f6e, 0x27ed86c, 0x1369671,
            0x36678a6, 0x1d6382f, 0x3b43d42, 0x1b913d9, 0x22e66d8, 0,
          },
        },
        {
          .x = {
            0, 0x3d5857a, 0x21033ac, 0x1bed27b, 0x2bcd8c8, 0x3dca518,
            0x394b264, 0x1d4c59a, 0x16efee4, 0x38261cb, 0x38e5426, 0,
          },
          .dt = {
            0, 0x39b0434, 0x4bda16, 0x225475e, 0x201ab85, 0x3f3dd8,
            0x1d3b648, 0x3195dab, 0x33b59ca, 0x1b2b3f1, 0x34ef35a, 0,
          },
          .y = {
            0, 0x554225, 0x38219a2, 0x219b8d3, 0x14b5ec4, 0x23c19a8,
            0xe9be87, 0x3fb6f05, 0x173773a, 0x3158979, 0x3f50a82, 0,
          },
        },
        {
          .x = {
            0, 0x164045e, 0x3b2c56a, 0x3c6ad21, 0x104aab6, 0x11496af,
            0x1e41c5d, 0x3572161, 0x3ff6f56, 0xca2f24, 0x1b6efc8, 0,
          },
          .dt = {
            0, 0x3b0e73, 0x2a38ad1, 0x230e852, 0x18e64f8, 0x3ce5194,
            0x33852f1, 0x2df29bb, 0x3209d88, 0x3b8cc15, 0xe68f3e, 0,
          },
          .y = {
            0, 0x3acb937, 0x2667273, 0x391b4c2, 0x3b6dede, 0x3ee7894,
            0x383e761, 0x34929d, 0xc42fe8, 0x26228fe, 0x35f6b5e, 0,
          },
        },
        {
          .x = {
            0, 0x17de973, 0x361143, 0x35e5b23, 0xe637b, 0x3fe0898,
            0xe618a6, 0x35cb223, 0x2bd339, 0xd7eb2e, 0x34d2423, 0,
          },
          .dt = {
            0, 0x3fbb7f1, 0xd36dd9, 0x3f77df5, 0x28418af, 0x5aab08,
            0x26ab4c9, 0x19f427d, 0xa60ff, 0x2481365, 0x1ed296f, 0,
          },
          .y = {
            0, 0xac37ab, 0x14f821e, 0x1109bc2, 0x172e885, 0x379b2c9,
            0x2def04b, 0xe76c84, 0x10de20f, 0x14b74d1, 0x125bd61, 0,
          },
        },
        {
          .x = {
            0, 0xe0d672, 0x2557fe1, 0x33428cb, 0x2a65ce6, 0x6b400f,
            0x2d9a97c, 0x1ffa163, 0x38fc4f9, 0x1374472, 0x43cea, 0,
          },
          .dt = {
            0, 0x2074317, 0x1ad72b2, 0x2010d5a, 0x1bf2cda, 0x2694f63,
            0x34ca7ac, 0x5c2c31, 0xf1ec8b, 0x25f574a, 0x1ad092, 0,
          },
          .y = {
            0, 0x2d0c780, 0xda3e29, 0x1aebe9a, 0x6b1b33, 0xfd4fae,
            0x127962f, 0x282b90e, 0x629049, 0x28e7d80, 0x14d1b78, 0,
          },
        },
        {
          .x = {
            0, 0x2648443, 0x23ba0d7, 0x31ea133, 0xc0142e, 0x32563d8,
            0x380cbca, 0x93a504, 0x1797fee, 0x15c3461, 0x1baed3c, 0,
          },
          .dt = {
            0, 0x2f1cf06, 0x216ca2, 0x3c42707, 0x14c153d, 0x1602e59,
            0x3426e97, 0x2104217, 0x3664f89, 0x2387c24, 0x3c276a5, 0,
          },
          .y = {
            0, 0x12aa53a, 0xd26204, 0x22acd56, 0x3ff45c8, 0x2c38501,
            0x21ff141, 0x3e6f963, 0x19da88c, 0x108a6a7, 0x22da2a5, 0,
          },
        },
        {
          .x = {
            0, 0x194d7cc, 0x33223ac, 0x2533f64, 0x219f8bb, 0x543263,
            0x236f648, 0xd9d222, 0x31e937b, 0x3c7583f, 0x1a1eaca, 0,
          },
          .dt = {
            0, 0x2d2c6c1, 0xd604ba, 0x264eca6, 0xd34607, 0x3588b9c,
            0x3a491cb, 0x1da1d8, 0x323d93b, 0x2c43651, 0x202a40a, 0,
          },
          .y = {
            0, 0x20c4cdd, 0x295a9cd, 0x15c8860, 0x22ef2d5, 0x222c240,
            0x2e4e2d1, 0xbea176, 0x365bf51, 0xf858ca, 0x294ed79, 0,
          },
        },
        {
          .x = {
            0, 0x14dd884, 0x12c9e33, 0x2d42122, 0x26f0b14, 0x1b9ea17,
            0x3779e94, 0x2562a88, 0xbe34f0, 0x192ead9, 0x89ec45, 0,
          },
          .dt = {
            0, 0x3ca3105, 0x30288c9, 0x3995bde, 0x795386, 0x1270116,
            0x3a8e6d2, 0x143f563, 0x15ad52, 0x38a2d02, 0x24daa5c, 0,
          },
          .y = {
            0, 0x1de5221, 0x172f820, 0x28c1b33, 0x8003c6, 0xe65926,
            0x188cd49, 0x3bb39fd, 0x1b9d8d7, 0x3d5020, 0x45742b, 0,
          },
        },
        {
          .x = {
            0, 0x1b40cbc, 0x1fc886e, 0x35f8918, 0x4ad538, 0x88cfee,
            0xdad340, 0x1b51fb0, 0x401665, 0xb82b1, 0x13af6c3, 0,
          },
          .dt = {
            0, 0x120880a, 0x15d2120, 0x8abfda, 0x110c710, 0xf21adf,
            0x1e72eef, 0x1285778, 0x8ceddc, 0x1559172, 0x2312acd, 0,
          },
          .y = {
            0, 0x1e3dcdb, 0x18ce3cb, 0x3b08774, 0xcff42, 0x3a41b4b,
            0x3dc7c05, 0x3da86bc, 0x83f2ce, 0x3d04d15, 0x31ade6f, 0,
          },
        },
        {
          .x = {
            0, 0x3bd6a2d, 0x3f51c16, 0x389ede, 0x3410d7e, 0xdba49,
            0x3b0c2dd, 0x37d1fed, 0x79084d, 0x3e00f0e, 0x1ff24ea, 0,
          },
          .dt = {
            0, 0x36580fc, 0x1c48fa9, 0x796c8d, 0x28176c8, 0x3dbe271,
            0x33a1f2, 0x1787d8c, 0xd5c744, 0xfb2dda, 0x2477293, 0,
          },
          .y = {
            0, 0x379e241, 0x1feb655, 0x362692e, 0x18a7163, 0x224ee73,
            0x18011e3, 0xc96c99, 0x386b934, 0x3ec531, 0x70840a, 0,
          },
        },
        {
          .x = {
            0, 0x42349a, 0x231743a, 0x34e2886, 0x367015e, 0x25a6fa0,
            0x1039803, 0x3afadc3, 0x28642e4, 0xb053a0, 0x384dd24, 0,
          },
          .dt = {
            0, 0x1c619a2, 0x8d14ad, 0x7b6876, 0x3f459ae, 0x182d60f,
            0x3e4af52, 0xe2730b, 0x10aa316, 0x1c426a, 0x13938a3, 0,
          },
          .y = {
            0, 0xaf4259, 0xd894f1, 0x152565, 0x1f4062, 0xcb9554,
            0x672d4d, 0x1925f1e, 0x935612, 0x202b0f2, 0x251a329, 0,
          },
        },
        {
          .x = {
            0, 0x28aaf07, 0x357061a, 0x39f926b, 0x1b814a2, 0x888c18,
            0x2d23fe4, 0x2fd9528, 0x2e8e707, 0x20069e1, 0x2a39f22, 0,
          },
          .dt = {
            0, 0x2d3d004, 0xa7a53e, 0x34b8a48, 0xa5b88b, 0x8dd67d,
            0xd394e4, 0x2b731b3, 0x73c0b8, 0xd992ad, 0x39fe6c6, 0,
          },
          .y = {
            0, 0x1b5263e, 0x3ee55da, 0x17c3c1f, 0x1952dbe, 0x2ff727e,
            0x5d31cb, 0x3557ac9, 0x74a614, 0x246210e, 0x303caeb, 0,
          },
        },
        {
          .x = {
            0, 0x3a3798b, 0x1182c2e, 0xef229e, 0x153db75, 0x3e23752,
            0x1b72274, 0x1b2d248, 0x39236ef, 0x35dff04, 0x261679a, 0,
          },
          .dt = {
            0, 0x3835012, 0x2b42c0e, 0xe41c7c, 0xdca92b, 0x386d306,
            0x197cd34, 0x1f75d0a, 0x22504e4, 0x2b1ab55, 0x14dfeb5, 0,
          },
          .y = {
            0, 0x57cd11, 0x2b809ae, 0x2f31ffc, 0x3f0ae37, 0x93c169,
            0x26131, 0x12dfdf8, 0x2b585d7, 0x5b2788, 0x270e65c, 0,
          },
        },
        {
          .x = {
            0, 0x315e4b6, 0xd7e8cc, 0x1462ee8, 0x1f7affe, 0x289ebf3,
            0x2daa6f8, 0x19dc553, 0x30672d9, 0x3db0c58, 0x3fd6e6d, 0,
          },
          .dt = {
            0, 0x1b92925, 0x267868e, 0x2bc9bef, 0x3809cc5, 0x31769a6,
            0x12ef439, 0x321136, 0x2cd81af, 0x1f74d8b, 0x22435c4, 0,
          },
          .y = {
            0, 0x11d34d7, 0x1064a61, 0x1e58bcc, 0x39fe27f, 0x2dbbae1,
            0x1a2ec16, 0x1993dfe, 0xe9dd4a, 0x3d07b24, 0x296d7f3, 0,
          },
        },
        {
          .x = {
            0, 0xfe2790, 0xe892d6, 0x15bd869, 0x124d184, 0x387f041,
            0x2912d87, 0x39eb14f, 0x2c2fcf9, 0xc31a89, 0x963208, 0,
          },
          .dt = {
            0, 0x1519046, 0x21a8288, 0x1e89af, 0x1d2c52d, 0x1cd5a2b,
            0x2fa34bb, 0x1771b6a, 0x147c7a5, 0x172f08d, 0x18cbd51, 0,
          },
          .y = {
            0, 0x96f2f7, 0x3a3f898, 0x2750e7c, 0xc7e76d, 0x33babe3,
            0x3c16c9f, 0x37aa3a2, 0x3032561, 0x3827e69, 0x9b5c2d, 0,
          },
        },
        {
          .x = {
            0, 0x20eef1a, 0x3c30e66, 0xd710f0, 0x248a6fa, 0x30c967f,
            0x3ce302c, 0xccd1f2, 0x197e993, 0x2ebaef3, 0xf2f019, 0,
          },
          .dt = {
            0, 0xe057f4, 0x1f50420, 0x135a47b, 0x225d8, 0x2014d14,
            0x2ae150f, 0x261e97c, 0x1932e2a, 0x880f54, 0xfdcba4, 0,
          },
          .y = {
            0, 0x3017cc0, 0x2a5110, 0x6d37e5, 0x283a64a, 0x1484b5,
            0x196f37b, 0x13de2d2, 0xda32d1, 0x392e0fc, 0x221d742, 0,
          },
        },
      },
    },
    {
      .table = {
        {
          .x = {
            0, 0x35fe901, 0x2743f3e, 0x3b80536, 0x21b386, 0x18f152d,
            0x196ba22, 0x11e9e36, 0x2e6f286, 0x36eb5d6, 0x1ef935f, 0,
          },
          .dt = {
            0, 0x11df0d1, 0x3e14615, 0x1a9aa54, 0x2875ff7, 0x292b2bb,
            0x3632d4d, 0x3a28ac3, 0x18261e, 0x4c2ddb, 0x37426ee, 0,
          },
          .y = {
            0, 0xac87cb, 0x208a864, 0x4f85a2, 0x2a5d88a, 0x3d6d1dd,
            0x2e71061, 0xca800a, 0x2407deb, 0x17abd58, 0xc62501, 0,
          },
        },
        {
          .x = {
            0, 0x3aed093, 0x127129, 0x2292952, 0x2a4d6dc, 0xcbd214,
            0x112f841, 0x1df8d6b, 0x1281b51, 0xd163e3, 0x3183603, 0,
          },
          .dt = {
            0, 0x6a6af5, 0x3c7596f, 0x2303d92, 0x7bb4de, 0x22e5bb0,
            0x3ccf944, 0x2e1d3f6, 0x1fbd8bc, 0x9851b4, 0x3233f62, 0,
          },
          .y = {
            0, 0x38b227, 0x1f6b7c, 0x72b604, 0x20e7c8e, 0x3bd6cc4,
            0xb86338, 0x607b90, 0xb47ed2, 0x4de572, 0x3fab39f, 0,
          },
        },
        {
          .x = {
            0, 0x7bcc46, 0x2fd7e67, 0x19cc479, 0x48d87e, 0x1e32cdc,
            0xb413f, 0x37cef98, 0x308ddc5, 0xf36b85, 0x10c9b00, 0,
          },
          .dt = {
            0, 0x1e8bca, 0x16a3b4a, 0x139b763, 0x318f6fa, 0x16ecf35,
            0x2eb9be6, 0x1c591eb, 0x1dfbb53, 0x25625d3, 0x399baaf, 0,
          },
          .y = {
            0, 0x2b6704b, 0x21c5085, 0x24b580f, 0x3fa439c, 0x2d70287,
            0x194a1c3, 0x30df043, 0x2d0a840, 0x1add393, 0x6d6688, 0,
          },
        },
        {
          .x = {
            0, 0x2e89bab, 0x890da0, 0x339e926, 0x3c5f59e, 0x1da780e,
            0x30b77ab, 0x93bfdb, 0x3507219, 0xddc262, 0x2e99c5, 0,
          },
          .dt = {
            0, 0x3790bf3, 0x15495ce, 0xe3efa1, 0x19b3423, 0xe1f075,
            0x1cffe32, 0xbf75d6, 0x2e1efc7, 0x15d8b7, 0x36aff6b, 0,
          },
          .y = {
            0, 0x3e80fb4, 0x72fbc2, 0x39502d6, 0x2d00e48, 0x227aedb,
            0x3410cd0, 0xb4eeb3, 0x18892a8, 0x3504f1b, 0x17cc6c4, 0,
          },
        },
        {
          .x = {
            0, 0xe60bff, 0x2c979d9, 0x1b1767, 0x3a0bcde, 0x2175d7d,
            0xc8e6d9, 0x167cb91, 0x2200d88, 0x1cf155a, 0x3e5e270, 0,
          },
          .dt = {
            0, 0x11212be, 0x1015507, 0x358501f, 0x13b5655, 0x3d52b60,
            0x9b8356, 0x2ed64d1, 0x396246d, 0x1e8c257, 0x18a4f4a, 0,
          },
          .y = {
            0, 0x188bf4, 0x32fe5f5, 0x9e87b1, 0x27eec6f, 0xf619fa,
            0x1a9d1b1, 0x36c77a1, 0x1116e84, 0x127bffd, 0xd4305b, 0,
          },
        },
        {
          .x = {
            0, 0x23fe686, 0x39b076f, 0x26372a0, 0x2e9a577, 0x29312bb,
            0x32a6d6c, 0x29f7b0, 0x13fb408, 0x1e79d5b, 0x1ff1d2c, 0,
          },
          .dt = {
            0, 0x113bd21, 0x7a08ed, 0xa8ea56, 0x3b62bcb, 0x3d7e3d1,
            0xe04b87, 0x10c0019, 0x17da9bc, 0x3502144, 0x228cfa7, 0,
          },
          .y = {
            0, 0x26546e0, 0x39fda01, 0x650c22, 0x3ac8c21, 0x132af4a,
            0x3bc5072, 0x3279814, 0x2c30e05, 0x2690852, 0x1e7e939, 0,
          },
        },
        {
          .x = {
            0, 0x3f4de3c, 0xee81b6, 0x31c45d3, 0x485207, 0x2d54c0b,
            0x22e44dc, 0x2d46130, 0x2fcee58, 0x1be150e, 0x202d8c4, 0,
          },
          .dt = {
            0, 0x18ecd22, 0x17ae926, 0x2f7f5b0, 0xcd8c4b, 0x393f150,
            0x31c5b08, 0x190ce1d, 0xdc72e4, 0x2625957, 0x2a8045b, 0,
          },
          .y = {
            0, 0x3885c57, 0x2b10f9b, 0xaf77c6, 0x1848a8c, 0x1dd2831,
            0x794231, 0x1e1e020, 0xdc8604, 0x1da9793, 0x3fb6328, 0,
          },
        },
        {
          .x = {
            0, 0x1d1cf29, 0x2e289d7, 0x1a83709, 0x2252d11, 0x3d6411c,
            0x3fd73ad, 0x2737d9c, 0x2ca9eba, 0x58f290, 0x3879a7c, 0,
          },
          .dt = {
            0, 0xc8ed55, 0x144a8cf, 0x20a0c, 0x1fd6b10, 0x689cd6,
            0x1d0eea, 0x1aadfdb, 0x75da9d, 0x13f3441, 0x37613dd, 0,
          },
          .y = {
            0, 0x357399d, 0x276752, 0xd5199f, 0x1bbd3a0, 0x39044f1,
            0xc5e83a, 0x1a99cdd, 0xdcb61f, 0x35b7272, 0x1184cff, 0,
          },
        },
        {
          .x = {
            0, 0x3a2435b, 0x2ee484b, 0x7cbb07, 0x18b3e31, 0x1cd6c22,
            0x1df9ba7, 0x10100b0, 0x1cb47fc, 0x980bc2, 0x38e7f39, 0,
          },
          .dt = {
            0, 0x39a3c84, 0x16bbf1b, 0x30639a6, 0x33a3a30, 0x2629ff5,
            0xccbe8c, 0x18cba1c, 0x1ce7453, 0x901515, 0x99efe1, 0,
          },
          .y = {
            0, 0x1ef39ce, 0x24e96f3, 0x3acf580, 0x3201d7e, 0x1d0de68,
            0x2cb19f3, 0x3505cc, 0xd96d33, 0x299dd6d, 0x1fa3ac7, 0,
          },
        },
        {
          .x = {
            0, 0x35e5061, 0x3b29489, 0x15a5454, 0x120682e, 0x27282d1,
            0x26ec442, 0x21ce37d, 0x2cf1802, 0x10fe952, 0xb2558c, 0,
          },
          .dt = {
            0, 0x3f2bef2, 0x3183f7b, 0x12a5cdf, 0x295adce, 0x12b738b,
            0xedfce, 0x1b62e, 0xc6f28b, 0x191d68f, 0x1e2a74d, 0,
          },
          .y = {
            0, 0x15536c8, 0x1c88676, 0x26a75b1, 0x35b75fb, 0x173eb88,
            0x2dd54c5, 0x162cbf1, 0x17fbd2, 0x1daef67, 0x15167dc, 0,
          },
        },
        {
          .x = {
            0, 0x22dcfea, 0x540134, 0x1caa84e, 0xba61aa, 0x1c70068,
            0x2e83a34, 0x15d4b1d, 0x21c6e92, 0x1a430ba, 0x210be9c, 0,
          },
          .dt = {
            0, 0xf9d61e, 0x3474cf3, 0x3d7a8b0, 0x250a2a4, 0x37046ad,
            0x23e45d1, 0x3593899, 0x663df7, 0x10a0e06, 0x1beb709, 0,
          },
          .y = {
            0, 0x2163326, 0x26bac62, 0x685312, 0x881acc, 0x183f4e2,
            0x3ec8a90, 0x3d030d1, 0x328ba50, 0x38fd0fc, 0xed5ed4, 0,
          },
        },
        {
          .x = {
            0, 0x1ebd890, 0x27161bf, 0x3279486, 0x2304a1e, 0x668ac3,
            0x327bc67, 0x3dea242, 0xaa7548, 0x2d360af, 0x3dac02c, 0,
          },
          .dt = {
            0, 0x276cf39, 0x51b201, 0x18130f9, 0x10622e0, 0x59b7a2,
            0xf3a49, 0x2add687, 0x21078fa, 0x314ee07, 0x3d2d76b, 0,
          },
          .y = {
            0, 0x211210c, 0x1048b49, 0xa35701, 0x7107b5, 0x1ecda2b,
            0x1b6a926, 0x3d439a7, 0x1a044ab, 0x1d983bd, 0x1407d5e, 0,
          },
        },
        {
          .x = {
            0, 0x243d680, 0x3beb552, 0x2891c48, 0x287efd3, 0x3b0f825,
            0x38caa11, 0x35dc81a, 0x3e0fe14, 0x632165, 0x768b30, 0,
          },
          .dt = {
            0, 0x273b8ff, 0x2b16c9f, 0x1fc2c93, 0x1a3f3ef, 0x1b138b9,
            0x370264f, 0x791a1e, 0x3aaba73, 0x40703d, 0x3a50f6b, 0,
          },
          .y = {
            0, 0xe86a2e, 0xe54b41, 0x3e34bcc, 0x3101bd1, 0x131caac,
            0x173d59e, 0x39c4c34, 0x481f61, 0xf4303a, 0x25a5bfd, 0,
          },
        },
        {
          .x = {
            0, 0x6166ae, 0x197a471, 0x289325a, 0x13959b4, 0x3c31d43,
            0xf2833a, 0x4fb591, 0x59f314, 0x731aff, 0x366b344, 0,
          },
          .dt = {
            0, 0x63af7f, 0x1fa5dd1, 0x2802e7d, 0xd22e2f, 0xa31528,
            0x3615725, 0x24090af, 0xfbcd4f, 0x1b88eaa, 0xfd71bc, 0,
          },
          .y = {
            0, 0x24bc0d8, 0x2ec1981, 0x2c5ca8d, 0x3edeaf7, 0x3fd0250,
            0x18320bf, 0x295ce3a, 0x32022e6, 0x22cd827, 0x1bc738a, 0,
          },
        },
        {
          .x = {
            0, 0x15f12f2, 0x337f0cd, 0x20aa879, 0x9fbeec, 0x1c28f22,
            0xbef02b, 0x253a408, 0x3238604, 0x18a8eab, 0x319f558, 0,
          },
          .dt = {
            0, 0x16a79f0, 0x152371b, 0x29ce7c4, 0x1002cdd, 0xa7755,
            0x1cf401e, 0x3f3f590, 0x16a3c63, 0x29fb307, 0x3fb1374, 0,
          },
          .y = {
            0, 0x3f20b18, 0xa7a2e1, 0x14cc854, 0x2613315, 0x1c05661,
            0x1e775f7, 0x224cb39, 0xe093dd, 0x35abd6c, 0x20de536, 0,
          },
        },
        {
          .x = {
            0, 0xe35d45, 0x38f90c, 0x283483, 0x1ee50a, 0x1e364f9,
            0x362414c, 0x156b1ed, 0x6fff6, 0x271f9ed, 0xffa45d, 0,
          },
          .dt = {
            0, 0x2f0c4c, 0x22ee1e1, 0x305a24f, 0x2b0d126, 0xa24164,
            0x13486cb, 0x820490, 0xffb00d, 0x1c414d8, 0x196878e, 0,
          },
          .y = {
            0, 0x156ae67, 0x27941ab, 0x19a3000, 0x3572ab5, 0x2b90ce3,
            0x136156c, 0x727496, 0xedae82, 0xfa5dfd, 0x16f293c, 0,
          },
        },
      },
    },
    {
      .table = {
        {
          .x = {
            0, 0x3321a95, 0x30b4a77, 0x33907a9, 0x1de8067, 0x841613,
            0x2f83d2, 0x1bd83e6, 0x1524483, 0x277922c, 0xcd77d0, 0,
          },
          .dt = {
            0, 0x1c6cbc6, 0x38f8eca, 0x1f8171a, 0x8a84e3, 0x3741a51,
            0x22190b4, 0x16cff5e, 0x2ae0fa0, 0x3e66c2f, 0xc13cfa, 0,
          },
          .y = {
            0, 0xcd2ebe, 0x21ac736, 0x2677442, 0x361f76, 0x132e1ae,
            0x333e7bd, 0x1da9fef, 0xafd793, 0x1f9c13a, 0x34b5916, 0,
          },
        },
        {
          .x = {
            0, 0xf3b2c0, 0x16cbccb, 0x1cffbb0, 0x1dba7a6, 0x2fc706d,
            0x145a053, 0x750c6c, 0x170c090, 0x319caa, 0x1fc4ceb, 0,
          },
          .dt = {
            0, 0x2d58554, 0x11b94a6, 0x72be6d, 0x2e06073, 0x14d8a1b,
            0x3c4088a, 0x34f27c8, 0x255e174, 0x6c8bad, 0x163a416, 0,
          },
          .y = {
            0, 0x160de28, 0x24202ce, 0x30931d4, 0xc11cd2, 0x1a13ee0,
            0xea6363, 0x3f3aee6, 0x263c832, 0x3330014, 0x12547ee, 0,
          },
        },
        {
          .x = {
            0, 0x884534, 0x58b8f7, 0x26e6659, 0x1b108d0, 0x18df97,
            0x251379b, 0x1000751, 0xcb0a03, 0x3c37fd2, 0x292f11f, 0,
          },
          .dt = {
            0, 0x34b7e79, 0x38713cb, 0xf7a976, 0x15044de, 0x3d727b,
            0x1f02470, 0x20f0f1f, 0xe0bae4, 0xee5c27, 0x3622d4f, 0,
          },
          .y = {
            0, 0x2f99c61, 0x38beb1b, 0x2b4260e, 0x3c51e81, 0x10358f7,
            0x15ed9fd, 0x9e54f3, 0x2577096, 0x5ff914, 0x763d19, 0,
          },
        },
        {
          .x = {
            0, 0x2dc03fc, 0x1f3944d, 0x1554033, 0xb75085, 0x28a773b,
            0x26a0431, 0xc53ef2, 0x29e0693, 0x20eb9d, 0x28a8186, 0,
          },
          .dt = {
            0, 0x1a42b75, 0x12a1c9a, 0x1dabaa8, 0xb7cbf6, 0x2e319af,
            0x2b4f86d, 0x1481b67, 0x260fa28, 0x129cef7, 0x20b4f15, 0,
          },
          .y = {
            0, 0x285bc2b, 0xe5aade, 0x356d640, 0x9fee7d, 0xddb16e,
            0x13e0729, 0x21abbb4, 0xaeb4eb, 0x8fd4cc, 0x1532ff6, 0,
          },
        },
        {
          .x = {
            0, 0x45d5d2, 0x2523f17, 0x1077568, 0x24ab8a5, 0x804b1b,
            0x8d65e2, 0x3f84910, 0x320ee64, 0x26b4d56, 0x11ff1cb, 0,
          },
          .dt = {
            0, 0x53fa51, 0x1213c48, 0x18305bd, 0x22d2bd0, 0x775408,
            0x3923281, 0x39a9a34, 0x3f1a960, 0x2e25d3b, 0x2c560ed, 0,
          },
          .y = {
            0, 0x251f03d, 0x39d6ba2, 0x35e2309, 0x1cbb0aa, 0xb7544a,
            0xffe530, 0x128ad0a, 0x17c7d40, 0x347668e, 0x3c399be, 0,
          },
        },
        {
          .x = {
            0, 0x129eb2e, 0x2420bd5, 0x356fd52, 0x184f0c1, 0xab660e,
            0xa282ef, 0x27276f3, 0x26de8e9, 0x2ce0b, 0x75f8d, 0,
          },
          .dt = {
            0, 0x33fcd14, 0x338388c, 0x146c5d4, 0x3515379, 0x1efac91,
            0x2f47998, 0x3917c92, 0x170ea71, 0x1184f7c, 0x875ed2, 0,
          },
          .y = {
            0, 0x3c587f4, 0x26d445, 0x2a567c3, 0x1be6358, 0x9e0909,
            0x1d65d2a, 0x3b8bd0f, 0x645e1d, 0x1436692, 0x33d0e67, 0,
          },
        },
        {
          .x = {
            0, 0x204588f, 0x2fb0b71, 0x16fee8a, 0x37604e8, 0x327d9ba,
            0x764080, 0x1e74286, 0x206ebad, 0x1601592, 0x3fda728, 0,
          },
          .dt = {
            0, 0x2432704, 0x19a148f, 0x154cc2a, 0x1cffd02, 0x39c3f01,
            0xc4aa9, 0x3309c6b, 0x6a58c3, 0x353751d, 0x1d8a4bc, 0,
          },
          .y = {
            0, 0x39e46ba, 0x1a4ca61, 0xba1507, 0x225fb90, 0x20ae696,
            0x2375799, 0xe28d6a, 0x6a73a0, 0xf2c02, 0x1290bcc, 0,
          },
        },
        {
          .x = {
            0, 0x1ea3c19, 0x81dc9e, 0x1a0b337, 0x1d7f3f4, 0x295a0aa,
            0x1ebff45, 0x956bf0, 0x17aae80, 0x5d8632, 0x3082c9a, 0,
          },
          .dt = {
            0, 0x2e156cd, 0x207b01, 0x3af4531, 0x360b9fc, 0x321ee9e,
            0x2bdf77, 0x13c6602, 0x374a146, 0x32744b7, 0x2723173, 0,
          },
          .y = {
            0, 0x22ad91f, 0x1ffcc65, 0x37b4f5c, 0x29c51ab, 0x3f9bd02,
            0x296aaf9, 0x2a58b82, 0x2c54e16, 0x2a7672c, 0x21486e2, 0,
          },
        },
        {
          .x = {
            0, 0x115406a, 0x3723cb2, 0x1eb8c9, 0x3ab5a08, 0x130c3cb,
            0x52d76, 0x24eb8d2, 0x385861f, 0xfdf882, 0x294e4d6, 0,
          },
          .dt = {
            0, 0x3a5daf5, 0x156cfa4, 0x3f62f39, 0x2de4863, 0x2d0dbfe,
            0x23952, 0x1211a9e, 0x388ff61, 0x373c759, 0x117ef53, 0,
          },
          .y = {
            0, 0x2a2f788, 0x2182dd0, 0x1090904, 0x36eb437, 0x3012d65,
            0x21c0716, 0x1437903, 0x23e0242, 0x3b9f9ae, 0xce8f2f, 0,
          },
        },
        {
          .x = {
            0, 0x13c3b88, 0x1b5537a, 0x28dbad2, 0xf1deee, 0x313e835,
            0x1c4a0d5, 0x83e575, 0x42aeae, 0x27c31a0, 0x1e14898, 0,
          },
          .dt = {
            0, 0x3ce3707, 0x335a7c2, 0x12c10ef, 0x28d7fed, 0x182b446,
            0x5fe9fc, 0x1ea6481, 0x2e23e91, 0x8c4f77, 0x33aefa4, 0,
          },
          .y = {
            0, 0x17899b, 0x29f9e15, 0x381698f, 0x481c59, 0x137de4b,
            0x25480e0, 0x1dd440b, 0x1308b65, 0x3493c23, 0x1134b8b, 0,
          },
        },
        {
          .x = {
            0, 0x301328, 0x25bc157, 0x161bcd8, 0x6b81f7, 0x27ec366,
            0x18f7abf, 0x1884c5b, 0x26c33cd, 0x2606c09, 0x2c69b7, 0,
          },
          .dt = {
            0, 0x15fb570, 0x250d45e, 0x5ef4b0, 0x3b7ddf9, 0x36d6260,
            0x168bfa4, 0x3201f67, 0x77dc0d, 0x2ca161, 0x124d0b9, 0,
          },
          .y = {
            0, 0x1ca8727, 0x294db9c, 0x2ce268e, 0x26dbe64, 0x364109b,
            0x2fb8996, 0x2274c9d, 0x3cddfd6, 0x2fe9408, 0x27d8bf4, 0,
          },
        },
        {
          .x = {
            0, 0x65a32c, 0x9e5ea3, 0xfa82b8, 0x33735dc, 0xd0decd,
            0x37e4fa0, 0x3f773c, 0xf41f87, 0x2d970b5, 0x1a66fbd, 0,
          },
          .dt = {
            0, 0x1cf5419, 0x2dac0e9, 0x294b93c, 0xa2015, 0x20c48b0,
            0x10454ee, 0xc6861d, 0x3734fff, 0x84d3bd, 0x34cb101, 0,
          },
          .y = {
            0, 0x3020dea, 0x1e7ad4e, 0x3c16e8d, 0xd8ab60, 0x90bc71,
            0x2d67dd2, 0x1f51415, 0x10534f7, 0x2856e58, 0x2940b35, 0,
          },
        },
        {
          .x = {
            0, 0x14bc9e6, 0x2214695, 0x1c54b00, 0x19312d1, 0x1ace064,
            0x615539, 0x165267f, 0x751f42, 0xd60433, 0x49571f, 0,
          },
          .dt = {
            0, 0x2727132, 0x334ddbe, 0x3035f3b, 0xdadd8, 0x3571f4c,
            0x22db60f, 0x349ebf7, 0x3351134, 0x9927e2, 0x280c2be, 0,
          },
          .y = {
            0, 0x5bac, 0xd1a346, 0x3a2f895, 0x3542ecd, 0x31a8fa6,
            0x33fa65d, 0x1693612, 0xff60a4, 0xc42002, 0x1bf14ab, 0,
          },
        },
        {
          .x = {
            0, 0x374e70c, 0x2a884e7, 0x35842dd, 0xf77b40, 0x348dc0f,
            0x1c152e3, 0x39fc297, 0x3e41076, 0xeffb13, 0x1298eeb, 0,
          },
          .dt = {
            0, 0xff7f63, 0x3fd2a, 0x1f011dc, 0x3fe59f0, 0x3b486b8,
            0x1a43473, 0xc1a580, 0x3700cf2, 0x23a7328, 0x2330049, 0,
          },
          .y = {
            0, 0x2c87159, 0x2e4eea7, 0x268d6e5, 0x500466, 0x3b61bb,
            0x164489e, 0x380b7c4, 0x11beb30, 0x33fe2bf, 0x1b0b9b2, 0,
          },
        },
        {
          .x = {
            0, 0x16c1f79, 0x176f5c6, 0x14edce6, 0x2002c64, 0x368117a,
            0xcc2c51, 0x31e05d2, 0x1ebed39, 0x1826838, 0x1505114, 0,
          },
          .dt = {
            0, 0x3acf28a, 0x23299fd, 0x659b4c, 0x1094ba1, 0x3cff25d,
            0x2d76be6, 0x312c122, 0x5710be, 0x2ffd7bb, 0xce2378, 0,
          },
          .y = {
            0, 0x39c33ab, 0x3e890f9, 0x418c3e, 0x27fadbb, 0x347cab0,
            0x23a0740, 0x1796bf4, 0xcd6e68, 0x3b8bc5, 0x2a54504, 0,
          },
        },
        {
          .x = {
            0, 0x37fcb1b, 0x16004b9, 0x1d18743, 0xbce648, 0xd78db6,
            0x35b1d65, 0x23bb620, 0x2fbc323, 0x1a9a586, 0x3b22577, 0,
          },
          .dt = {
            0, 0xfdfd2c, 0x8245dd, 0x3645fdd, 0x25ebad0, 0x30eb1db,
            0x229b719, 0x1ddb6d1, 0xa8798e, 0x39a21da, 0x36f57b4, 0,
          },
          .y = {
            0, 0x82fb15, 0x3487d6, 0x3d1c2c9, 0x2c9e7ad, 0x187be10,
            0x2e9b6ba, 0x15b8f89, 0x243ae4c, 0x328bb11, 0xb12a9, 0,
          },
        },
      },
    },
    {
      .table = {
        {
          .x = {
            0, 0x1663380, 0x9c3a78, 0x1ba0201, 0x2d90546, 0x1c5afe5,
            0x33c57ab, 0x22daff6, 0x2a9dfbe, 0x28ef998, 0x60d74d, 0,
          },
          .dt = {
            0, 0x368d6e0, 0x27e3ec1, 0x1edc55c, 0x174b367, 0x177d3e3,
            0x3fe3f95, 0x1a1528c, 0x3d0b364, 0x36f45d8, 0x27e4fb5, 0,
          },
          .y = {
            0, 0x33f2824, 0x2cb48c0, 0x1b32c9b, 0x3595cf8, 0x39c074b,
            0x2ff8009, 0x1a6311, 0x1b5ee70, 0xdf2f29, 0xf1a06d, 0,
          },
        },
        {
          .x = {
            0, 0x4de2bb, 0x2746b84, 0x2a8b258, 0xec80dd, 0x24b8036,
            0x3cb3935, 0x364ae61, 0x114d8ad, 0x189aacd, 0x34244a9, 0,
          },
          .dt = {
            0, 0x18b0d28, 0x3bb116d, 0x6e82fc, 0x3fc747a, 0x304de82,
            0x1a72c8c, 0x1d30fe5, 0x2fbf1d5, 0xcdd270, 0x205ef9c, 0,
          },
          .y = {
            0, 0x18f4a78, 0x2d0cd3e, 0x1199e07, 0x1363597, 0x18f467,
            0x1dfa043, 0x10f3e0c, 0xc4bfe1, 0x29b87b4, 0x33a7be7, 0,
          },
        },
        {
          .x = {
            0, 0xb369c5, 0x1632733, 0x2e8680f, 0x21cfa28, 0x381141b,
            0x376c23a, 0x1273078, 0x4369f1, 0x575c5, 0x2f8dc16, 0,
          },
          .dt = {
            0, 0x1ad0ce2, 0x3124a52, 0x2e3a34, 0x3ed7343, 0x27eef57,
            0x1a43f2c, 0x21193db, 0x3fd904, 0xab9776, 0x1d42152, 0,
          },
          .y = {
            0, 0xf4a9ae, 0x37e55e5, 0x272059e, 0x24a79b6, 0x376c8a,
            0x1cb6510, 0x323f093, 0x170f65e, 0x32b0d96, 0x1a5c2c6, 0,
          },
        },
        {
          .x = {
            0, 0x1ce1b, 0xc01eca, 0x279eadf, 0x37ceb0c, 0xb7c8d4,
            0x3ad2f0, 0x344137e, 0x3c3556f, 0x224f31d, 0x2a4f634, 0,
          },
          .dt = {
            0, 0x230e553, 0x292fd80, 0x364e744, 0x15d5b6, 0x205f8a0,
            0xf24073, 0x19471f2, 0x1b90b97, 0x1828ef3, 0x3033be9, 0,
          },
          .y = {
            0, 0x30c03ed, 0x3838efa, 0x2c0e01f, 0x3882638, 0x2058a85,
            0x11bd1da, 0x3b94f34, 0x162fe44, 0x3cee7aa, 0x1517aec, 0,
          },
        },
        {
          .x = {
            0, 0x2f6a040, 0x623b2, 0x1fdf9a7, 0x712345, 0x3a0e99a,
            0x17e2821, 0x1a6e7b5, 0x1b204f2, 0x106d7b2, 0x2b79a60, 0,
          },
          .dt = {
            0, 0x32045c3, 0x1f1e2ed, 0x1c8964, 0x9448ed, 0xba2d18,
            0xfd01d0, 0xc75330, 0x2d7cc61, 0x1dd5084, 0x389a335, 0,
          },
          .y = {
            0, 0x696d91, 0x3f4e8a8, 0x1efdc62, 0x39359e3, 0x39a70db,
            0x293bd2e, 0x2e4c889, 0x3e3017f, 0xbd4cf3, 0xf63b4b, 0,
          },
        },
        {
          .x = {
            0, 0x16e30ab, 0x3b66579, 0x357c449, 0x402047, 0x1595280,
            0x991bd4, 0x2f5cf7d, 0x314bb15, 0x19fa0c7, 0x2750b8d, 0,
          },
          .dt = {
            0, 0x28cb55, 0xea4373, 0x2ec6b, 0xae2fcf, 0xd2bf6e,
            0x3aa79f6, 0x3d14684, 0x2464b65, 0x2bf3b47, 0xa6e369, 0,
          },
          .y = {
            0, 0x213bf35, 0x37b19d1, 0x15095a3, 0x31e0da8, 0x146569c,
            0x20c44b8, 0x31d6ccc, 0x3ba0c57, 0x24162af, 0x28c4ee, 0,
          },
        },
        {
          .x = {
            0, 0x2e18721, 0x3351332, 0x1707880, 0x2ae8750, 0x1b2bba1,
            0x57dd62, 0x38d0141, 0xe97911, 0x273b133, 0x34db3b5, 0,
          },
          .dt = {
            0, 0x49c6d8, 0x4c65e4, 0x1d64149, 0x6a87d5, 0x15afd59,
            0x2064a93, 0xf80380, 0x2dd4c6d, 0x25167c, 0x2491a08, 0,
          },
          .y = {
            0, 0x288d18a, 0x2e61da6, 0x250a878, 0x11d3870, 0x85fd68,
            0x1c180f8, 0x2fcd1de, 0x3313c03, 0x32455e, 0x109787a, 0,
          },
        },
        {
          .x = {
            0, 0x6b9c9d, 0x3d00674, 0x10a73fc, 0x30fda83, 0x139185c,
            0x43e082, 0x3c67915, 0x208192a, 0x25e451, 0x258a566, 0,
          },
          .dt = {
            0, 0x1eced50, 0x2fba908, 0x24cb5f1, 0x26e71ed, 0x35f1a16,
            0x757df0, 0x3314110, 0x27265c0, 0x10e7282, 0x1b2fb0e, 0,
          },
          .y = {
            0, 0x3d2a04f, 0x1314c36, 0x131c7a3, 0x1882ef3, 0x1a0a5e8,
            0x1919356, 0xa5616a, 0x1eea31d, 0x2c216b3, 0x18ba4aa, 0,
          },
        },
        {
          .x = {
            0, 0x2eecff0, 0x181971f, 0x12edf62, 0x3d6dca5, 0x2f1a768,
            0x1b4b84f, 0xeb8ff2, 0x1ee870a, 0x1daf31d, 0x1929218, 0,
          },
          .dt = {
            0, 0x22138fa, 0x3743418, 0x3754221, 0x1d2d615, 0x32cd2b9,
            0xe6efb0, 0x3965d7a, 0xdb06b6, 0x26c29a, 0x2b66056, 0,
          },
          .y = {
            0, 0x2968e29, 0x2410dd9, 0x4aa8ef, 0x2940073, 0x2f93c25,
            0x19d20cf, 0x1603915, 0x25e4e0f, 0x2b063c, 0x2a7390e, 0,
          },
        },
        {
          .x = {
            0, 0xb57819, 0x1804f41, 0x32bb607, 0x327ff64, 0x17adb8c,
            0x789a23, 0x353c021, 0x196b696, 0x2dfb7cc, 0x2dd8c7e, 0,
          },
          .dt = {
            0, 0x37a2bdd, 0x392bb41, 0x34790ed, 0x2eb5b12, 0x13943a4,
            0x15a5f40, 0xb831b5, 0x30ea094, 0x2a8e922, 0x218a382, 0,
          },
          .y = {
            0, 0x734ead, 0x1e0aa23, 0x25f9c21, 0x2653674, 0x277f303,
            0x32167aa, 0x316cbf3, 0x25b38f0, 0x2ef80b7, 0x3cfb960, 0,
          },
        },
        {
          .x = {
            0, 0x27b4a46, 0x147ad41, 0x17d59ac, 0xafe100, 0x3704e6d,
            0x2b8be90, 0x17fe9ee, 0x22de9e5, 0x1b6cfd0, 0x3c2e60d, 0,
          },
          .dt = {
            0, 0x3fe6936, 0x35a6545, 0x369fc71, 0x2d2df4e, 0x2858fb6,
            0x295b66, 0x12a21d9, 0x2e5a512, 0x1a0b5c4, 0x36e9cf6, 0,
          },
          .y = {
            0, 0x31e059e, 0x339e4d1, 0xc92e5b, 0x1ab881f, 0x357985f,
            0x3f8c14a, 0x27fe7cb, 0x26d021a, 0x52b6a0, 0x18a446b, 0,
          },
        },
        {
          .x = {
            0, 0x2eab0dd, 0x1fdbcf1, 0x26be8d5, 0x248b48f, 0x3c28eb5,
            0xb057ba, 0x1c81a86, 0x1eac18, 0xf697b8, 0x334b103, 0,
          },
          .dt = {
            0, 0x2a1427f, 0x12543cc, 0x23e64e, 0xcadedf, 0x2581bdc,
            0x1c6dee1, 0x10edb96, 0x1c131c9, 0x896be8, 0x342ccce, 0,
          },
          .y = {
            0, 0x1d8e833, 0x3ae9f2a, 0x3bfdcf, 0x3e03b8b, 0x2b6c6a4,
            0xca8e7a, 0xbd92ba, 0x2edc80c, 0x2586dbf, 0x164527c, 0,
          },
        },
        {
          .x = {
            0, 0x2d2205c, 0xb5bacd, 0x6ea15c, 0xe1eb83, 0x12784fe,
            0x2d776cc, 0x2d33c43, 0x1ab5b1f, 0x320cdee, 0x3fd72ca, 0,
          },
          .dt = {
            0, 0x23165f3, 0x48f672, 0x37c8b58, 0x464c74, 0xb819ea,
            0x1a5a289, 0x1ef821e, 0x49d352, 0xdf5312, 0x67ceb8, 0,
          },
          .y = {
            0, 0xeb7070, 0x1bcd642, 0x3e50ee0, 0x336a937, 0x39477f4,
            0x16919db, 0xf7d1a4, 0xb5e7ad, 0x3809791, 0x2a75d4c, 0,
          },
        },
        {
          .x = {
            0, 0x1de0cbf, 0x4c6c92, 0x27a97d8, 0x1c1e22f, 0x30e7677,
            0xe27480, 0x1a5dcf1, 0x23c6b40, 0x74d032, 0x16093d4, 0,
          },
          .dt = {
            0, 0x3dc37b8, 0x2d89e94, 0x38d16c9, 0x3572753, 0x3037644,
            0x348bdb4, 0xb5990, 0x3dead9a, 0x1648fe3, 0x3967bd5, 0,
          },
          .y = {
            0, 0x1acfabe, 0x21911e6, 0x367e61d, 0x2eb8fab, 0x37c7660,
            0x36ae5b1, 0x437a8c, 0x248ad46, 0x3efe641, 0x3cb605b, 0,
          },
        },
        {
          .x = {
            0, 0x63a1c8, 0x285c9ff, 0xc82f49, 0x146bdd0, 0x8e2bc7,
            0x3241e22, 0x35e4839, 0x144ebbc, 0x8322b7, 0x1c1178a, 0,
          },
          .dt = {
            0, 0x1c23a80, 0x192db94, 0x1a18cc6, 0x2eb7015, 0x35e83fa,
            0x183ced8, 0x2e259c3, 0x2db1b88, 0x312c4ce, 0x2588509, 0,
          },
          .y = {
            0, 0x1fe096f, 0x1500328, 0x3d97e46, 0x30fe95, 0x3f640cd,
            0x3b836e9, 0x230b3e3, 0x1a8f125, 0x3e23de6, 0x2e66752, 0,
          },
        },
        {
          .x = {
            0, 0x3e79b25, 0x2ca71b7, 0x2b2ea3c, 0xde7ac4, 0x3026d10,
            0x2bce79e, 0x1153866, 0x3e5a80, 0x22b9a37, 0x3e9c59, 0,
          },
          .dt = {
            0, 0x10020e2, 0x3423da6, 0xe0197, 0x1a2eb25, 0x263506e,
            0x26d121b, 0x2a62e35, 0x57494b, 0x8f4aca, 0x1848d1, 0,
          },
          .y = {
            0, 0x20100d6, 0x2330974, 0x3402585, 0x172cfd6, 0x275a21c,
            0x213e87c, 0x29989f2, 0x155e437, 0x96a378, 0x3a674eb, 0,
          },
        },
      },
    },
  },
};

void compute_comb_set(
  sabs_comb_set_t *result, const affine_pt_narrow_t *base_pt) {
  sabs_comb_set_wide_t result_t;
  // Temporarily store the individual teeth for computing the individual comb
  // entries.
  teeth_set_t teeth_sets[COMB_COUNT];

  // For signed all bits set representation, we start with the all bits positive
  // representation. We build those points incrementally in this array.
  extended_pt_wide_t everything_pts[COMB_COUNT];

  // Used to hold the current point during the doubling operation. We start the
  // iteration at 1 because we do special things at the start of a comb, and the
  // base point needs to have those done differently.
  projective_pt_wide_t temp;
  extended_pt_wide_t temp_ext;
  // Start at 2^1 * b
  affine_double_extended(&temp_ext, base_pt);
  // We'll need the projective point for the start of the loop.
  extended_to_projective_wide(&temp, &temp_ext);
  // We start the teeth sets out negative because we're going to alternate signs
  // during the gray code walk. Note that we don't put the base point into the
  // tooth set. (Because we never need it as a delta to walk SABS gray code)
  extended_to_readd_wide_neg(&teeth_sets[0].teeth[0], &temp_ext);
  // Initialize the first everything point
  affine_narrow_to_extended(&everything_pts[0], base_pt);
  for (int i = 1; i < COMB_TEETH * COMB_COUNT; ++i) {
    // We stop 2 before: 1 So we can get an extended result, and a second
    // because SABS gives us the highest bit for free.
    for (int j = 0; j < COMB_SEPARATION - 2; ++j) {
      projective_double(&temp, &temp);
    }
    projective_double_extended(&temp_ext, &temp);
    if (i % COMB_TEETH == 0) {
      // If this is the first tooth in the comb, just copy it into the
      // everything array.
      copy_extended_pt_wide(&everything_pts[i/COMB_TEETH], &temp_ext);
    } else {
      // Otherwise add it to the already in progress everything point.
      extended_add_extended(
        &everything_pts[i/COMB_TEETH], &everything_pts[i/COMB_TEETH],
        &temp_ext);
    }
    // If this isn't the very last iteration
    if (i != COMB_TEETH * COMB_COUNT - 1) {
      if (i % COMB_TEETH != COMB_TEETH - 1) {
        // If it isn't the last tooth of the comb
        extended_double_extended(&temp_ext, &temp_ext);
        extended_to_readd_wide_neg(
          &teeth_sets[i/COMB_TEETH].teeth[i % COMB_TEETH], &temp_ext);
        // We'll need the projective point for the loop.
        extended_to_projective_wide(&temp, &temp_ext);
      } else {
        // Double to match the double we would have done if this weren't the
        // last tooth.
        extended_to_projective_wide(&temp, &temp_ext);
        projective_double(&temp, &temp);
      }
    }
  }
  // We now have all the precomputation necessary to walk the gray codes and
  // compute the table entries.
  int entry = COMB_TABLE_SIZE - 1;
  // The everything point goes in the very last slot in the comb table
  for (int i = 0; i < COMB_COUNT; ++i) {
    extended_to_projective_wide(
      &result_t.combs[i].table[entry], &everything_pts[i]);
  }
  for (int i = 1; i < COMB_TABLE_SIZE; ++i) {
    int j;
    // This makes entry walk a gray code. The gray code starts at
    // COMB_TABLE_SIZE - 1 rather than starting at 0. At the end of the loop, j
    // is the bit that was toggled.
    for (j = 0; j < COMB_TEETH - 1; ++j) {
      int bit = 1 << (j + 1);
      int mask = bit - 1;
      int half_bit = bit >> 1;
      if ((i & mask) == half_bit) {
        entry ^= half_bit;
        break;
      }
    }
    // Toggle bit j in every comb.
    for (int k = 0; k < COMB_COUNT; ++k) {
      extended_readd_wide_extended(
        &everything_pts[k], &everything_pts[k], &teeth_sets[k].teeth[j]);
      // To prepare for future toggling, negate the tooth.
      negate_extended_pt_readd_wide(
        &teeth_sets[k].teeth[j], &teeth_sets[k].teeth[j]);
      // Store the result while also going from extended to projective
      extended_to_projective_wide(
        &result_t.combs[k].table[entry], &everything_pts[k]);
    }
  }
  // Finally reduce to get the ultimate result
  reduce_comb_set(result, &result_t);
}

// Note not const. Stomps on source
// Leaves source with invalid z-values, but if they are all set to one, it will
// be correct again. Shrinks combs down to affine with a single inversion. Also
// narrows and reduces to 10 limbs.
void reduce_comb_set(sabs_comb_set_t *result, sabs_comb_set_wide_t *source) {
  residue_wide_t z_left;
  residue_wide_t z_right;

  residue_wide_t z_inv;

  // This is an interview question that is actually applicable. For every entry
  // in the table we want to multiply it by every z entry but its own. We
  // do it in two passes, accumulating z values as we go from left to right and
  // right to left.
  copy_wide(&z_left, &source->combs[0].table[0].z);
  for (int i = 1; i < COMB_TABLE_SIZE * COMB_COUNT; ++i) {
    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &z_left);
    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &z_left);
    mul_wide(&z_left, &z_left,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].z);
  }

  copy_wide(&z_right,
      &source->combs[COMB_COUNT - 1].table[COMB_TABLE_SIZE - 1].z);

  for (int i = COMB_TABLE_SIZE * COMB_COUNT - 2; i >= 0; --i) {
    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &z_right);
    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &z_right);
    mul_wide(&z_right, &z_right,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].z);
  }

  // Now that all the entries have common denominators, perform a single
  // inversion
  invert_wide(&z_inv, &z_right);
  for (int i = 0; i < COMB_TABLE_SIZE * COMB_COUNT; ++i) {
    residue_wide_t xy;

    // Divide by the common z and narrow.
    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &z_inv);
    narrow(&result->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
           &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x);

    mul_wide(&source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
             &z_inv);
    narrow(&result->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y,
           &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y);

    // The combs store DT rather than T for faster additions.
    mul_wide(&xy,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].x,
             &source->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].y);
    mul_wide_const(&xy, &xy, D);
    narrow(
      &result->combs[i / COMB_TABLE_SIZE].table[i % COMB_TABLE_SIZE].dt,
      &xy);
  }
}

void scalar_comb_multiply(
  projective_pt_wide_t *result, const sabs_comb_set_t * __restrict comb,
  const scalar_t * __restrict n) {

  scalar_t sabs_n;
  convert_to_sabs(&sabs_n, n);

  extended_pt_wide_t temp;
  extended_affine_pt_readd_narrow_t table_pt;

  // Start with the highest bits because we double the accumulator
  for (int i = COMB_SEPARATION - 1; i >= 0; --i) {
    if (i != COMB_SEPARATION - 1) {
      extended_double_extended(&temp, &temp);
    }
    for (int j = 0; j < COMB_COUNT; ++j) {
      int entry = 0;

      // extract the specific bits for this comb entry
      for (int k = 0; k < COMB_TEETH; ++k) {
        int bit = i + COMB_SEPARATION * (k + COMB_TEETH * j);
        if (bit < SCALAR_BITS) {
          entry |= ((sabs_n.limbs[bit / SCALAR_LIMB_BITS] >>
              (bit % SCALAR_LIMB_BITS)) & 1) << k;
        }
      }

      // The highest bit is the sign bit.
      int32_t invert = (entry >> (COMB_TEETH - 1)) - 1;
      entry ^= invert;

      constant_time_extended_affine_narrow_lookup(
        &table_pt, entry & COMB_LOOKUP_MASK, COMB_TABLE_SIZE,
        comb->combs[j].table);

      constant_time_cond_extended_affine_negate(&table_pt, invert);

      if (i == (COMB_SEPARATION - 1) && j == 0) {
        affine_readd_to_extended(&temp, &table_pt);
      } else {
        extended_readd_affine_narrow_extended(
          &temp, &temp, &table_pt);
      }
    }
  }

  extended_to_projective_wide(result, &temp);
  explicit_bzero(&sabs_n, sizeof(sabs_n));
  explicit_bzero(&table_pt, sizeof(table_pt));
  explicit_bzero(&temp, sizeof(temp));
}

void scalar_comb_multiply_unsafe(
  projective_pt_wide_t *result, const sabs_comb_set_t * __restrict comb,
  const scalar_t * __restrict n) {

  scalar_t sabs_n;
  convert_to_sabs(&sabs_n, n);

  extended_pt_wide_t temp;
  extended_affine_pt_readd_narrow_t table_pt;

  for (int i = COMB_SEPARATION - 1; i >= 0; --i) {
    if (i != COMB_SEPARATION - 1) {
      extended_double_extended(&temp, &temp);
    }
    for (int j = 0; j < COMB_COUNT; ++j) {
      int entry = 0;

      for (int k = 0; k < COMB_TEETH; ++k) {
        int bit = i + COMB_SEPARATION * (k + COMB_TEETH * j);
        if (bit < SCALAR_BITS) {
          entry |= ((sabs_n.limbs[bit / SCALAR_LIMB_BITS] >>
              (bit % SCALAR_LIMB_BITS)) & 1) << k;
        }
      }

      int32_t invert = (entry >> (COMB_TEETH - 1)) - 1;
      entry ^= invert;

      copy_extended_affine_pt_readd_narrow(
        &table_pt, &comb->combs[j].table[entry & COMB_LOOKUP_MASK]);

      if (invert) {
        negate_extended_affine_pt_readd_narrow(&table_pt, &table_pt);
      }

      if (i == (COMB_SEPARATION - 1) && j == 0) {
        affine_readd_to_extended(&temp, &table_pt);
      } else {
        extended_readd_affine_narrow_extended(
          &temp, &temp, &table_pt);
      }
    }
  }

  extended_to_projective_wide(result, &temp);
  explicit_bzero(&sabs_n, sizeof(sabs_n));
  explicit_bzero(&table_pt, sizeof(table_pt));
  explicit_bzero(&temp, sizeof(temp));
}
