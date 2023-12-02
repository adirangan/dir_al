/* This holds a mem-list item */
struct mitem
{
  struct mitem *parent;
  struct mitem *child;
  void *item;
  size_t size;
};

/* This holds a list of mitems */
struct mlist
{
  struct mitem *first;
  struct mitem *last;
  int length;
  size_t size;
};

/* This holds a link-listed item */
struct litem
{
  struct litem *parent;
  struct litem *child;
  void *item;
};

/* This holds a link-list of litems */
struct llist
{
  struct litem *first;
  struct litem *last;
  int length;
};

/* This holds a llitem */
struct llitem
{
  struct llitem *parent;
  struct llitem *kidl;
  struct llitem *kidr;
  int nright;
  int nleft;
  void *item;
};

/* This holds a histogram */
struct hist
{
  int nbins;
  double max;
  double min;
  double *data;
  double datasum;
  double ndatum;
};

/* Here are RNUM variables -- readonable values are int RNEXT=2;int RPREV=1;int RADD=3;int RPOW=16; */
static long long int POW2RPOWPLUSRADD=65539LL;
static long long int POW22RPOWMINUSONE=2147483648LL;

/* Here are nodes and weights for 32 bit legendre quadrature */
static int GLOBAL_QUADRATURE_N=32;
static double GLOBAL_QUADRATURE_W[32]={0.70186100094700966004E-02, 0.16274394730905670605E-01, 0.25392065309262059456E-01, 0.34273862913021433103E-01, 0.42835898022226680657E-01, 0.50998059262376176196E-01, 0.58684093478535547145E-01, 0.65822222776361846838E-01, 0.72345794108848506225E-01, 0.78193895787070306472E-01, 0.83311924226946755222E-01, 0.87652093004403811143E-01, 0.91173878695763884713E-01, 0.93844399080804565639E-01, 0.95638720079274859419E-01, 0.96540088514727800567E-01, 0.96540088514727800567E-01, 0.95638720079274859419E-01, 0.93844399080804565639E-01, 0.91173878695763884713E-01, 0.87652093004403811143E-01, 0.83311924226946755222E-01, 0.78193895787070306472E-01, 0.72345794108848506225E-01, 0.65822222776361846838E-01, 0.58684093478535547145E-01, 0.50998059262376176196E-01, 0.42835898022226680657E-01, 0.34273862913021433103E-01, 0.25392065309262059456E-01, 0.16274394730905670605E-01, 0.70186100094700966004E-02};
static double GLOBAL_QUADRATURE_Z[32]={-0.99726386184948156354E+00, -0.98561151154526833540E+00, -0.96476225558750643077E+00, -0.93490607593773968917E+00, -0.89632115576605212397E+00, -0.84936761373256997013E+00, -0.79448379596794240696E+00, -0.73218211874028968039E+00, -0.66304426693021520098E+00, -0.58771575724076232904E+00, -0.50689990893222939002E+00, -0.42135127613063534536E+00, -0.33186860228212764978E+00, -0.23928736225213707454E+00, -0.14447196158279649349E+00, -0.48307665687738316235E-01, 0.48307665687738316235E-01, 0.14447196158279649349E+00, 0.23928736225213707454E+00, 0.33186860228212764978E+00, 0.42135127613063534536E+00, 0.50689990893222939002E+00, 0.58771575724076232904E+00, 0.66304426693021520098E+00, 0.73218211874028968039E+00, 0.79448379596794240696E+00, 0.84936761373256997013E+00, 0.89632115576605212397E+00, 0.93490607593773968917E+00, 0.96476225558750643077E+00, 0.98561151154526833540E+00, 0.99726386184948156354E+00};

/* Here are preambles for standard and 512 color figure generation */
static char FIG_PREAMBLE[128]="#FIG 3.2\nLandscape\nCenter\nInches\nLetter\n100.00\nSingle\n-2\n1200 2\n";
static char FIG_PREAMBLE_COLOR_0[8192]="0 32 #000000\n0 33 #000000\n0 34 #010101\n0 35 #010101\n0 36 #020202\n0 37 #020202\n0 38 #030303\n0 39 #030303\n0 40 #040404\n0 41 #040404\n0 42 #050505\n0 43 #050505\n0 44 #060606\n0 45 #060606\n0 46 #070707\n0 47 #070707\n0 48 #080808\n0 49 #080808\n0 50 #090909\n0 51 #090909\n0 52 #0a0a0a\n0 53 #0a0a0a\n0 54 #0b0b0b\n0 55 #0b0b0b\n0 56 #0c0c0c\n0 57 #0c0c0c\n0 58 #0d0d0d\n0 59 #0d0d0d\n0 60 #0e0e0e\n0 61 #0e0e0e\n0 62 #0f0f0f\n0 63 #0f0f0f\n0 64 #101010\n0 65 #101010\n0 66 #111111\n0 67 #111111\n0 68 #121212\n0 69 #121212\n0 70 #131313\n0 71 #131313\n0 72 #141414\n0 73 #141414\n0 74 #151515\n0 75 #151515\n0 76 #161616\n0 77 #161616\n0 78 #171717\n0 79 #171717\n0 80 #181818\n0 81 #181818\n0 82 #191919\n0 83 #191919\n0 84 #1a1a1a\n0 85 #1a1a1a\n0 86 #1b1b1b\n0 87 #1b1b1b\n0 88 #1c1c1c\n0 89 #1c1c1c\n0 90 #1d1d1d\n0 91 #1d1d1d\n0 92 #1e1e1e\n0 93 #1e1e1e\n0 94 #1f1f1f\n0 95 #1f1f1f\n0 96 #202020\n0 97 #202020\n0 98 #212121\n0 99 #212121\n0 100 #222222\n0 101 #222222\n0 102 #232323\n0 103 #232323\n0 104 #242424\n0 105 #242424\n0 106 #252525\n0 107 #252525\n0 108 #262626\n0 109 #262626\n0 110 #272727\n0 111 #272727\n0 112 #282828\n0 113 #282828\n0 114 #292929\n0 115 #292929\n0 116 #2a2a2a\n0 117 #2a2a2a\n0 118 #2b2b2b\n0 119 #2b2b2b\n0 120 #2c2c2c\n0 121 #2c2c2c\n0 122 #2d2d2d\n0 123 #2d2d2d\n0 124 #2e2e2e\n0 125 #2e2e2e\n0 126 #2f2f2f\n0 127 #2f2f2f\n0 128 #303030\n0 129 #303030\n0 130 #313131\n0 131 #313131\n0 132 #323232\n0 133 #323232\n0 134 #333333\n0 135 #333333\n0 136 #343434\n0 137 #343434\n0 138 #353535\n0 139 #353535\n0 140 #363636\n0 141 #363636\n0 142 #373737\n0 143 #373737\n0 144 #383838\n0 145 #383838\n0 146 #393939\n0 147 #393939\n0 148 #3a3a3a\n0 149 #3a3a3a\n0 150 #3b3b3b\n0 151 #3b3b3b\n0 152 #3c3c3c\n0 153 #3c3c3c\n0 154 #3d3d3d\n0 155 #3d3d3d\n0 156 #3e3e3e\n0 157 #3e3e3e\n0 158 #3f3f3f\n0 159 #3f3f3f\n0 160 #404040\n0 161 #404040\n0 162 #414141\n0 163 #414141\n0 164 #424242\n0 165 #424242\n0 166 #434343\n0 167 #434343\n0 168 #444444\n0 169 #444444\n0 170 #454545\n0 171 #454545\n0 172 #464646\n0 173 #464646\n0 174 #474747\n0 175 #474747\n0 176 #484848\n0 177 #484848\n0 178 #494949\n0 179 #494949\n0 180 #4a4a4a\n0 181 #4a4a4a\n0 182 #4b4b4b\n0 183 #4b4b4b\n0 184 #4c4c4c\n0 185 #4c4c4c\n0 186 #4d4d4d\n0 187 #4d4d4d\n0 188 #4e4e4e\n0 189 #4e4e4e\n0 190 #4f4f4f\n0 191 #4f4f4f\n0 192 #505050\n0 193 #505050\n0 194 #515151\n0 195 #515151\n0 196 #525252\n0 197 #525252\n0 198 #535353\n0 199 #535353\n0 200 #545454\n0 201 #545454\n0 202 #555555\n0 203 #555555\n0 204 #565656\n0 205 #565656\n0 206 #575757\n0 207 #575757\n0 208 #585858\n0 209 #585858\n0 210 #595959\n0 211 #595959\n0 212 #5a5a5a\n0 213 #5a5a5a\n0 214 #5b5b5b\n0 215 #5b5b5b\n0 216 #5c5c5c\n0 217 #5c5c5c\n0 218 #5d5d5d\n0 219 #5d5d5d\n0 220 #5e5e5e\n0 221 #5e5e5e\n0 222 #5f5f5f\n0 223 #5f5f5f\n0 224 #606060\n0 225 #606060\n0 226 #616161\n0 227 #616161\n0 228 #626262\n0 229 #626262\n0 230 #636363\n0 231 #636363\n0 232 #646464\n0 233 #646464\n0 234 #656565\n0 235 #656565\n0 236 #666666\n0 237 #666666\n0 238 #676767\n0 239 #676767\n0 240 #686868\n0 241 #686868\n0 242 #696969\n0 243 #696969\n0 244 #6a6a6a\n0 245 #6a6a6a\n0 246 #6b6b6b\n0 247 #6b6b6b\n0 248 #6c6c6c\n0 249 #6c6c6c\n0 250 #6d6d6d\n0 251 #6d6d6d\n0 252 #6e6e6e\n0 253 #6e6e6e\n0 254 #6f6f6f\n0 255 #6f6f6f\n0 256 #707070\n0 257 #707070\n0 258 #717171\n0 259 #717171\n0 260 #727272\n0 261 #727272\n0 262 #737373\n0 263 #737373\n0 264 #747474\n0 265 #747474\n0 266 #757575\n0 267 #757575\n0 268 #767676\n0 269 #767676\n0 270 #777777\n0 271 #777777\n0 272 #787878\n0 273 #787878\n0 274 #797979\n0 275 #797979\n0 276 #7a7a7a\n0 277 #7a7a7a\n0 278 #7b7b7b\n0 279 #7b7b7b\n0 280 #7c7c7c\n0 281 #7c7c7c\n0 282 #7d7d7d\n0 283 #7d7d7d\n0 284 #7e7e7e\n0 285 #7e7e7e\n0 286 #7f7f7f\n0 287 #7f7f7f\n0 288 #808080\n0 289 #808080\n0 290 #818181\n0 291 #818181\n0 292 #828282\n0 293 #828282\n0 294 #838383\n0 295 #838383\n0 296 #848484\n0 297 #848484\n0 298 #858585\n0 299 #858585\n0 300 #868686\n0 301 #868686\n0 302 #878787\n0 303 #878787\n0 304 #888888\n0 305 #888888\n0 306 #898989\n0 307 #898989\n0 308 #8a8a8a\n0 309 #8a8a8a\n0 310 #8b8b8b\n0 311 #8b8b8b\n0 312 #8c8c8c\n0 313 #8c8c8c\n0 314 #8d8d8d\n0 315 #8d8d8d\n0 316 #8e8e8e\n0 317 #8e8e8e\n0 318 #8f8f8f\n0 319 #8f8f8f\n0 320 #909090\n0 321 #909090\n0 322 #919191\n0 323 #919191\n0 324 #929292\n0 325 #929292\n0 326 #939393\n0 327 #939393\n0 328 #949494\n0 329 #949494\n0 330 #959595\n0 331 #959595\n0 332 #969696\n0 333 #969696\n0 334 #979797\n0 335 #979797\n0 336 #989898\n0 337 #989898\n0 338 #999999\n0 339 #999999\n0 340 #9a9a9a\n0 341 #9a9a9a\n0 342 #9b9b9b\n0 343 #9b9b9b\n0 344 #9c9c9c\n0 345 #9c9c9c\n0 346 #9d9d9d\n0 347 #9d9d9d\n0 348 #9e9e9e\n0 349 #9e9e9e\n0 350 #9f9f9f\n0 351 #9f9f9f\n0 352 #a0a0a0\n0 353 #a0a0a0\n0 354 #a1a1a1\n0 355 #a1a1a1\n0 356 #a2a2a2\n0 357 #a2a2a2\n0 358 #a3a3a3\n0 359 #a3a3a3\n0 360 #a4a4a4\n0 361 #a4a4a4\n0 362 #a5a5a5\n0 363 #a5a5a5\n0 364 #a6a6a6\n0 365 #a6a6a6\n0 366 #a7a7a7\n0 367 #a7a7a7\n0 368 #a8a8a8\n0 369 #a8a8a8\n0 370 #a9a9a9\n0 371 #a9a9a9\n0 372 #aaaaaa\n0 373 #aaaaaa\n0 374 #ababab\n0 375 #ababab\n0 376 #acacac\n0 377 #acacac\n0 378 #adadad\n0 379 #adadad\n0 380 #aeaeae\n0 381 #aeaeae\n0 382 #afafaf\n0 383 #afafaf\n0 384 #b0b0b0\n0 385 #b0b0b0\n0 386 #b1b1b1\n0 387 #b1b1b1\n0 388 #b2b2b2\n0 389 #b2b2b2\n0 390 #b3b3b3\n0 391 #b3b3b3\n0 392 #b4b4b4\n0 393 #b4b4b4\n0 394 #b5b5b5\n0 395 #b5b5b5\n0 396 #b6b6b6\n0 397 #b6b6b6\n0 398 #b7b7b7\n0 399 #b7b7b7\n0 400 #b8b8b8\n0 401 #b8b8b8\n0 402 #b9b9b9\n0 403 #b9b9b9\n0 404 #bababa\n0 405 #bababa\n0 406 #bbbbbb\n0 407 #bbbbbb\n0 408 #bcbcbc\n0 409 #bcbcbc\n0 410 #bdbdbd\n0 411 #bdbdbd\n0 412 #bebebe\n0 413 #bebebe\n0 414 #bfbfbf\n0 415 #bfbfbf\n0 416 #c0c0c0\n0 417 #c0c0c0\n0 418 #c1c1c1\n0 419 #c1c1c1\n0 420 #c2c2c2\n0 421 #c2c2c2\n0 422 #c3c3c3\n0 423 #c3c3c3\n0 424 #c4c4c4\n0 425 #c4c4c4\n0 426 #c5c5c5\n0 427 #c5c5c5\n0 428 #c6c6c6\n0 429 #c6c6c6\n0 430 #c7c7c7\n0 431 #c7c7c7\n0 432 #c8c8c8\n0 433 #c8c8c8\n0 434 #c9c9c9\n0 435 #c9c9c9\n0 436 #cacaca\n0 437 #cacaca\n0 438 #cbcbcb\n0 439 #cbcbcb\n0 440 #cccccc\n0 441 #cccccc\n0 442 #cdcdcd\n0 443 #cdcdcd\n0 444 #cecece\n0 445 #cecece\n0 446 #cfcfcf\n0 447 #cfcfcf\n0 448 #d0d0d0\n0 449 #d0d0d0\n0 450 #d1d1d1\n0 451 #d1d1d1\n0 452 #d2d2d2\n0 453 #d2d2d2\n0 454 #d3d3d3\n0 455 #d3d3d3\n0 456 #d4d4d4\n0 457 #d4d4d4\n0 458 #d5d5d5\n0 459 #d5d5d5\n0 460 #d6d6d6\n0 461 #d6d6d6\n0 462 #d7d7d7\n0 463 #d7d7d7\n0 464 #d8d8d8\n0 465 #d8d8d8\n0 466 #d9d9d9\n0 467 #d9d9d9\n0 468 #dadada\n0 469 #dadada\n0 470 #dbdbdb\n0 471 #dbdbdb\n0 472 #dcdcdc\n0 473 #dcdcdc\n0 474 #dddddd\n0 475 #dddddd\n0 476 #dedede\n0 477 #dedede\n0 478 #dfdfdf\n0 479 #dfdfdf\n0 480 #e0e0e0\n0 481 #e0e0e0\n0 482 #e1e1e1\n0 483 #e1e1e1\n0 484 #e2e2e2\n0 485 #e2e2e2\n0 486 #e3e3e3\n0 487 #e3e3e3\n0 488 #e4e4e4\n0 489 #e4e4e4\n0 490 #e5e5e5\n0 491 #e5e5e5\n0 492 #e6e6e6\n0 493 #e6e6e6\n0 494 #e7e7e7\n0 495 #e7e7e7\n0 496 #e8e8e8\n0 497 #e8e8e8\n0 498 #e9e9e9\n0 499 #e9e9e9\n0 500 #eaeaea\n0 501 #eaeaea\n0 502 #ebebeb\n0 503 #ebebeb\n0 504 #ececec\n0 505 #ececec\n0 506 #ededed\n0 507 #ededed\n0 508 #eeeeee\n0 509 #eeeeee\n0 510 #efefef\n0 511 #efefef\n0 512 #f0f0f0\n0 513 #f0f0f0\n0 514 #f1f1f1\n0 515 #f1f1f1\n0 516 #f2f2f2\n0 517 #f2f2f2\n0 518 #f3f3f3\n0 519 #f3f3f3\n0 520 #f4f4f4\n0 521 #f4f4f4\n0 522 #f5f5f5\n0 523 #f5f5f5\n0 524 #f6f6f6\n0 525 #f6f6f6\n0 526 #f7f7f7\n0 527 #f7f7f7\n0 528 #f8f8f8\n0 529 #f8f8f8\n0 530 #f9f9f9\n0 531 #f9f9f9\n0 532 #fafafa\n0 533 #fafafa\n0 534 #fbfbfb\n0 535 #fbfbfb\n0 536 #fcfcfc\n0 537 #fcfcfc\n0 538 #fdfdfd\n0 539 #fdfdfd\n0 540 #fefefe\n0 541 #fefefe\n0 542 #ffffff\n0 543 #ffffff\n";
static char FIG_PREAMBLE_COLOR_5[8192]="0 32 #0000ff\n0 33 #0002ff\n0 34 #0004ff\n0 35 #0006ff\n0 36 #0008ff\n0 37 #000aff\n0 38 #000cff\n0 39 #000eff\n0 40 #000fff\n0 41 #0011ff\n0 42 #0013ff\n0 43 #0015ff\n0 44 #0017fe\n0 45 #0019fe\n0 46 #001bfe\n0 47 #001dfe\n0 48 #001ffe\n0 49 #0021fd\n0 50 #0023fd\n0 51 #0025fd\n0 52 #0027fc\n0 53 #0029fc\n0 54 #002bfc\n0 55 #002dfb\n0 56 #002ffb\n0 57 #0031fb\n0 58 #0033fa\n0 59 #0034fa\n0 60 #0036fa\n0 61 #0038f9\n0 62 #003af9\n0 63 #003cf8\n0 64 #003ef8\n0 65 #0040f7\n0 66 #0041f7\n0 67 #0043f6\n0 68 #0045f6\n0 69 #0047f5\n0 70 #0048f5\n0 71 #004af4\n0 72 #004cf4\n0 73 #004ef3\n0 74 #004ff3\n0 75 #0051f2\n0 76 #0053f2\n0 77 #0055f1\n0 78 #0056f0\n0 79 #0058f0\n0 80 #005aef\n0 81 #005bef\n0 82 #005dee\n0 83 #005eed\n0 84 #0060ed\n0 85 #0062ec\n0 86 #0063eb\n0 87 #0065eb\n0 88 #0066ea\n0 89 #0068e9\n0 90 #0069e9\n0 91 #006be8\n0 92 #006ce7\n0 93 #006ee7\n0 94 #006fe6\n0 95 #0071e5\n0 96 #0072e4\n0 97 #0074e4\n0 98 #0075e3\n0 99 #0076e2\n0 100 #0078e1\n0 101 #0079e1\n0 102 #007be0\n0 103 #007cdf\n0 104 #007ddf\n0 105 #007fde\n0 106 #0080dd\n0 107 #0081dc\n0 108 #0082dc\n0 109 #0084db\n0 110 #0085da\n0 111 #0086d9\n0 112 #0087d8\n0 113 #0089d8\n0 114 #008ad7\n0 115 #008bd6\n0 116 #008cd5\n0 117 #008dd5\n0 118 #008ed4\n0 119 #0090d3\n0 120 #0091d2\n0 121 #0092d2\n0 122 #0093d1\n0 123 #0094d0\n0 124 #0095cf\n0 125 #0096ce\n0 126 #0097ce\n0 127 #0098cd\n0 128 #0099cc\n0 129 #009acb\n0 130 #009bcb\n0 131 #009cca\n0 132 #009dc9\n0 133 #009ec8\n0 134 #009fc8\n0 135 #00a0c7\n0 136 #00a1c6\n0 137 #00a2c5\n0 138 #00a3c5\n0 139 #00a4c4\n0 140 #00a5c3\n0 141 #00a6c2\n0 142 #00a7c1\n0 143 #00a7c1\n0 144 #00a8c0\n0 145 #00a9bf\n0 146 #00aabf\n0 147 #00abbe\n0 148 #00acbd\n0 149 #00acbc\n0 150 #00adbc\n0 151 #00aebb\n0 152 #00afba\n0 153 #00b0b9\n0 154 #00b0b9\n0 155 #00b1b8\n0 156 #00b2b7\n0 157 #00b3b6\n0 158 #00b3b6\n0 159 #00b4b5\n0 160 #00b5b4\n0 161 #00b5b4\n0 162 #00b6b3\n0 163 #00b7b2\n0 164 #00b8b1\n0 165 #00b8b1\n0 166 #00b9b0\n0 167 #00baaf\n0 168 #00baae\n0 169 #00bbae\n0 170 #00bcad\n0 171 #00bdac\n0 172 #00bdab\n0 173 #00beaa\n0 174 #00bfaa\n0 175 #00c0a9\n0 176 #00c0a8\n0 177 #00c1a7\n0 178 #00c2a6\n0 179 #00c3a5\n0 180 #00c3a4\n0 181 #00c4a3\n0 182 #00c5a3\n0 183 #00c6a2\n0 184 #00c6a1\n0 185 #00c7a0\n0 186 #00c89f\n0 187 #00c99e\n0 188 #00c99d\n0 189 #00ca9c\n0 190 #00cb9b\n0 191 #00cc9a\n0 192 #00cd99\n0 193 #00cd98\n0 194 #00ce97\n0 195 #00cf96\n0 196 #00d095\n0 197 #00d093\n0 198 #00d192\n0 199 #00d291\n0 200 #00d390\n0 201 #00d38f\n0 202 #00d48e\n0 203 #00d58d\n0 204 #00d68c\n0 205 #00d78a\n0 206 #00d789\n0 207 #00d888\n0 208 #00d987\n0 209 #00da86\n0 210 #00da84\n0 211 #00db83\n0 212 #00dc82\n0 213 #00dd80\n0 214 #00dd7f\n0 215 #00de7e\n0 216 #00df7d\n0 217 #00e07b\n0 218 #00e07a\n0 219 #00e178\n0 220 #00e277\n0 221 #00e376\n0 222 #00e374\n0 223 #00e473\n0 224 #00e571\n0 225 #00e570\n0 226 #00e66f\n0 227 #00e76d\n0 228 #00e86c\n0 229 #00e86a\n0 230 #00e969\n0 231 #00ea67\n0 232 #00ea66\n0 233 #00eb64\n0 234 #00ec62\n0 235 #00ec61\n0 236 #00ed5f\n0 237 #00ee5e\n0 238 #00ee5c\n0 239 #00ef5a\n0 240 #00ef59\n0 241 #00f057\n0 242 #00f155\n0 243 #00f154\n0 244 #00f252\n0 245 #00f250\n0 246 #00f34f\n0 247 #00f44d\n0 248 #00f44b\n0 249 #00f549\n0 250 #00f548\n0 251 #00f646\n0 252 #00f644\n0 253 #00f742\n0 254 #00f740\n0 255 #00f83f\n0 256 #00f83d\n0 257 #00f93b\n0 258 #00f939\n0 259 #00f937\n0 260 #00fa35\n0 261 #00fa33\n0 262 #00fb32\n0 263 #00fb30\n0 264 #00fb2e\n0 265 #00fc2c\n0 266 #00fc2a\n0 267 #00fc28\n0 268 #00fd26\n0 269 #00fd24\n0 270 #00fd22\n0 271 #00fd20\n0 272 #00fe1e\n0 273 #00fe1c\n0 274 #00fe1a\n0 275 #00fe18\n0 276 #00fe16\n0 277 #00ff14\n0 278 #00ff12\n0 279 #00ff10\n0 280 #00ff0f\n0 281 #00ff0d\n0 282 #00ff0b\n0 283 #00ff09\n0 284 #00ff07\n0 285 #00ff05\n0 286 #00ff03\n0 287 #00ff01\n0 288 #01ff00\n0 289 #03ff00\n0 290 #05ff00\n0 291 #07ff00\n0 292 #09ff00\n0 293 #0bff00\n0 294 #0dff00\n0 295 #0fff00\n0 296 #10ff00\n0 297 #12ff00\n0 298 #14ff00\n0 299 #16fe00\n0 300 #18fe00\n0 301 #1afe00\n0 302 #1cfe00\n0 303 #1efe00\n0 304 #20fd00\n0 305 #22fd00\n0 306 #24fd00\n0 307 #26fd00\n0 308 #28fc00\n0 309 #2afc00\n0 310 #2cfc00\n0 311 #2efb00\n0 312 #30fb00\n0 313 #32fb00\n0 314 #33fa00\n0 315 #35fa00\n0 316 #37f900\n0 317 #39f900\n0 318 #3bf900\n0 319 #3df800\n0 320 #3ff800\n0 321 #40f700\n0 322 #42f700\n0 323 #44f600\n0 324 #46f600\n0 325 #48f500\n0 326 #49f500\n0 327 #4bf400\n0 328 #4df400\n0 329 #4ff300\n0 330 #50f200\n0 331 #52f200\n0 332 #54f100\n0 333 #55f100\n0 334 #57f000\n0 335 #59ef00\n0 336 #5aef00\n0 337 #5cee00\n0 338 #5eee00\n0 339 #5fed00\n0 340 #61ec00\n0 341 #62ec00\n0 342 #64eb00\n0 343 #66ea00\n0 344 #67ea00\n0 345 #69e900\n0 346 #6ae800\n0 347 #6ce800\n0 348 #6de700\n0 349 #6fe600\n0 350 #70e500\n0 351 #71e500\n0 352 #73e400\n0 353 #74e300\n0 354 #76e300\n0 355 #77e200\n0 356 #78e100\n0 357 #7ae000\n0 358 #7be000\n0 359 #7ddf00\n0 360 #7ede00\n0 361 #7fdd00\n0 362 #80dd00\n0 363 #82dc00\n0 364 #83db00\n0 365 #84da00\n0 366 #86da00\n0 367 #87d900\n0 368 #88d800\n0 369 #89d700\n0 370 #8ad700\n0 371 #8cd600\n0 372 #8dd500\n0 373 #8ed400\n0 374 #8fd300\n0 375 #90d300\n0 376 #91d200\n0 377 #92d100\n0 378 #93d000\n0 379 #95d000\n0 380 #96cf00\n0 381 #97ce00\n0 382 #98cd00\n0 383 #99cd00\n0 384 #9acc00\n0 385 #9bcb00\n0 386 #9cca00\n0 387 #9dc900\n0 388 #9ec900\n0 389 #9fc800\n0 390 #a0c700\n0 391 #a1c600\n0 392 #a2c600\n0 393 #a3c500\n0 394 #a3c400\n0 395 #a4c300\n0 396 #a5c300\n0 397 #a6c200\n0 398 #a7c100\n0 399 #a8c000\n0 400 #a9c000\n0 401 #aabf00\n0 402 #aabe00\n0 403 #abbd00\n0 404 #acbd00\n0 405 #adbc00\n0 406 #aebb00\n0 407 #aeba00\n0 408 #afba00\n0 409 #b0b900\n0 410 #b1b800\n0 411 #b1b800\n0 412 #b2b700\n0 413 #b3b600\n0 414 #b4b500\n0 415 #b4b500\n0 416 #b5b400\n0 417 #b6b300\n0 418 #b6b300\n0 419 #b7b200\n0 420 #b8b100\n0 421 #b9b000\n0 422 #b9b000\n0 423 #baaf00\n0 424 #bbae00\n0 425 #bcad00\n0 426 #bcac00\n0 427 #bdac00\n0 428 #beab00\n0 429 #bfaa00\n0 430 #bfa900\n0 431 #c0a800\n0 432 #c1a700\n0 433 #c1a700\n0 434 #c2a600\n0 435 #c3a500\n0 436 #c4a400\n0 437 #c5a300\n0 438 #c5a200\n0 439 #c6a100\n0 440 #c7a000\n0 441 #c89f00\n0 442 #c89e00\n0 443 #c99d00\n0 444 #ca9c00\n0 445 #cb9b00\n0 446 #cb9a00\n0 447 #cc9900\n0 448 #cd9800\n0 449 #ce9700\n0 450 #ce9600\n0 451 #cf9500\n0 452 #d09400\n0 453 #d19300\n0 454 #d29200\n0 455 #d29100\n0 456 #d39000\n0 457 #d48e00\n0 458 #d58d00\n0 459 #d58c00\n0 460 #d68b00\n0 461 #d78a00\n0 462 #d88900\n0 463 #d88700\n0 464 #d98600\n0 465 #da8500\n0 466 #db8400\n0 467 #dc8200\n0 468 #dc8100\n0 469 #dd8000\n0 470 #de7f00\n0 471 #df7d00\n0 472 #df7c00\n0 473 #e07b00\n0 474 #e17900\n0 475 #e17800\n0 476 #e27600\n0 477 #e37500\n0 478 #e47400\n0 479 #e47200\n0 480 #e57100\n0 481 #e66f00\n0 482 #e76e00\n0 483 #e76c00\n0 484 #e86b00\n0 485 #e96900\n0 486 #e96800\n0 487 #ea6600\n0 488 #eb6500\n0 489 #eb6300\n0 490 #ec6200\n0 491 #ed6000\n0 492 #ed5e00\n0 493 #ee5d00\n0 494 #ef5b00\n0 495 #ef5a00\n0 496 #f05800\n0 497 #f05600\n0 498 #f15500\n0 499 #f25300\n0 500 #f25100\n0 501 #f34f00\n0 502 #f34e00\n0 503 #f44c00\n0 504 #f44a00\n0 505 #f54800\n0 506 #f54700\n0 507 #f64500\n0 508 #f64300\n0 509 #f74100\n0 510 #f74000\n0 511 #f83e00\n0 512 #f83c00\n0 513 #f93a00\n0 514 #f93800\n0 515 #fa3600\n0 516 #fa3400\n0 517 #fa3300\n0 518 #fb3100\n0 519 #fb2f00\n0 520 #fb2d00\n0 521 #fc2b00\n0 522 #fc2900\n0 523 #fc2700\n0 524 #fd2500\n0 525 #fd2300\n0 526 #fd2100\n0 527 #fe1f00\n0 528 #fe1d00\n0 529 #fe1b00\n0 530 #fe1900\n0 531 #fe1700\n0 532 #ff1500\n0 533 #ff1300\n0 534 #ff1100\n0 535 #ff0f00\n0 536 #ff0e00\n0 537 #ff0c00\n0 538 #ff0a00\n0 539 #ff0800\n0 540 #ff0600\n0 541 #ff0400\n0 542 #ff0200\n0 543 #ff0000\n";
static char FIG_PREAMBLE_COLOR_7[8192]="0 32 #0000ff\n0 33 #0002ff\n0 34 #0005ff\n0 35 #0007ff\n0 36 #000aff\n0 37 #000cff\n0 38 #000fff\n0 39 #0011ff\n0 40 #0014ff\n0 41 #0016ff\n0 42 #0019ff\n0 43 #001bff\n0 44 #001eff\n0 45 #0020ff\n0 46 #0023ff\n0 47 #0025ff\n0 48 #0028ff\n0 49 #002aff\n0 50 #002dff\n0 51 #002fff\n0 52 #0032ff\n0 53 #0034ff\n0 54 #0037ff\n0 55 #0039ff\n0 56 #003cff\n0 57 #003eff\n0 58 #0041ff\n0 59 #0043ff\n0 60 #0046ff\n0 61 #0048ff\n0 62 #004bff\n0 63 #004dff\n0 64 #0050ff\n0 65 #0052ff\n0 66 #0055ff\n0 67 #0057ff\n0 68 #005aff\n0 69 #005cff\n0 70 #005fff\n0 71 #0061ff\n0 72 #0064ff\n0 73 #0066ff\n0 74 #0069ff\n0 75 #006bff\n0 76 #006eff\n0 77 #0070ff\n0 78 #0073ff\n0 79 #0075ff\n0 80 #0078ff\n0 81 #007aff\n0 82 #007dff\n0 83 #007fff\n0 84 #0082ff\n0 85 #0084ff\n0 86 #0087ff\n0 87 #0089ff\n0 88 #008cff\n0 89 #008eff\n0 90 #0091ff\n0 91 #0093ff\n0 92 #0096ff\n0 93 #0098ff\n0 94 #009bff\n0 95 #009dff\n0 96 #00a0ff\n0 97 #00a2ff\n0 98 #00a5ff\n0 99 #00a7ff\n0 100 #00aaff\n0 101 #00acff\n0 102 #00afff\n0 103 #00b1ff\n0 104 #00b4ff\n0 105 #00b6ff\n0 106 #00b9ff\n0 107 #00bbff\n0 108 #00beff\n0 109 #00c0ff\n0 110 #00c3ff\n0 111 #00c5ff\n0 112 #00c8ff\n0 113 #00caff\n0 114 #00cdff\n0 115 #00cfff\n0 116 #00d2ff\n0 117 #00d4ff\n0 118 #00d7ff\n0 119 #00d9ff\n0 120 #00dcff\n0 121 #00deff\n0 122 #00e1ff\n0 123 #00e3ff\n0 124 #00e6ff\n0 125 #00e8ff\n0 126 #00ebff\n0 127 #00edff\n0 128 #00f0ff\n0 129 #00f2ff\n0 130 #00f5ff\n0 131 #00f7ff\n0 132 #00faff\n0 133 #00fcff\n0 134 #00ffff\n0 135 #00fffd\n0 136 #00fffb\n0 137 #00fff8\n0 138 #00fff6\n0 139 #00fff3\n0 140 #00fff1\n0 141 #00ffee\n0 142 #00ffec\n0 143 #00ffe9\n0 144 #00ffe7\n0 145 #00ffe4\n0 146 #00ffe2\n0 147 #00ffdf\n0 148 #00ffdd\n0 149 #00ffda\n0 150 #00ffd8\n0 151 #00ffd5\n0 152 #00ffd3\n0 153 #00ffd0\n0 154 #00ffce\n0 155 #00ffcb\n0 156 #00ffc9\n0 157 #00ffc6\n0 158 #00ffc4\n0 159 #00ffc1\n0 160 #00ffbf\n0 161 #00ffbc\n0 162 #00ffba\n0 163 #00ffb7\n0 164 #00ffb5\n0 165 #00ffb2\n0 166 #00ffb0\n0 167 #00ffad\n0 168 #00ffab\n0 169 #00ffa8\n0 170 #00ffa6\n0 171 #00ffa3\n0 172 #00ffa1\n0 173 #00ff9e\n0 174 #00ff9c\n0 175 #00ff99\n0 176 #00ff97\n0 177 #00ff94\n0 178 #00ff92\n0 179 #00ff8f\n0 180 #00ff8d\n0 181 #00ff8a\n0 182 #00ff88\n0 183 #00ff85\n0 184 #00ff83\n0 185 #00ff80\n0 186 #00ff7e\n0 187 #00ff7b\n0 188 #00ff79\n0 189 #00ff76\n0 190 #00ff74\n0 191 #00ff71\n0 192 #00ff6f\n0 193 #00ff6c\n0 194 #00ff6a\n0 195 #00ff67\n0 196 #00ff65\n0 197 #00ff62\n0 198 #00ff60\n0 199 #00ff5d\n0 200 #00ff5b\n0 201 #00ff58\n0 202 #00ff56\n0 203 #00ff53\n0 204 #00ff51\n0 205 #00ff4e\n0 206 #00ff4c\n0 207 #00ff49\n0 208 #00ff47\n0 209 #00ff44\n0 210 #00ff42\n0 211 #00ff3f\n0 212 #00ff3d\n0 213 #00ff3a\n0 214 #00ff38\n0 215 #00ff35\n0 216 #00ff33\n0 217 #00ff30\n0 218 #00ff2e\n0 219 #00ff2b\n0 220 #00ff29\n0 221 #00ff26\n0 222 #00ff24\n0 223 #00ff21\n0 224 #00ff1f\n0 225 #00ff1c\n0 226 #00ff1a\n0 227 #00ff17\n0 228 #00ff15\n0 229 #00ff12\n0 230 #00ff10\n0 231 #00ff0d\n0 232 #00ff0b\n0 233 #00ff08\n0 234 #00ff06\n0 235 #00ff03\n0 236 #00ff01\n0 237 #01ff00\n0 238 #04ff00\n0 239 #06ff00\n0 240 #09ff00\n0 241 #0bff00\n0 242 #0eff00\n0 243 #10ff00\n0 244 #13ff00\n0 245 #15ff00\n0 246 #18ff00\n0 247 #1aff00\n0 248 #1dff00\n0 249 #1fff00\n0 250 #22ff00\n0 251 #24ff00\n0 252 #27ff00\n0 253 #29ff00\n0 254 #2cff00\n0 255 #2eff00\n0 256 #31ff00\n0 257 #33ff00\n0 258 #36ff00\n0 259 #38ff00\n0 260 #3bff00\n0 261 #3dff00\n0 262 #40ff00\n0 263 #42ff00\n0 264 #45ff00\n0 265 #47ff00\n0 266 #4aff00\n0 267 #4cff00\n0 268 #4fff00\n0 269 #51ff00\n0 270 #54ff00\n0 271 #56ff00\n0 272 #59ff00\n0 273 #5bff00\n0 274 #5eff00\n0 275 #60ff00\n0 276 #63ff00\n0 277 #65ff00\n0 278 #68ff00\n0 279 #6aff00\n0 280 #6dff00\n0 281 #6fff00\n0 282 #72ff00\n0 283 #74ff00\n0 284 #77ff00\n0 285 #79ff00\n0 286 #7cff00\n0 287 #7eff00\n0 288 #81ff00\n0 289 #83ff00\n0 290 #86ff00\n0 291 #88ff00\n0 292 #8bff00\n0 293 #8dff00\n0 294 #90ff00\n0 295 #92ff00\n0 296 #95ff00\n0 297 #97ff00\n0 298 #9aff00\n0 299 #9cff00\n0 300 #9fff00\n0 301 #a1ff00\n0 302 #a4ff00\n0 303 #a6ff00\n0 304 #a9ff00\n0 305 #abff00\n0 306 #aeff00\n0 307 #b0ff00\n0 308 #b3ff00\n0 309 #b5ff00\n0 310 #b8ff00\n0 311 #baff00\n0 312 #bdff00\n0 313 #bfff00\n0 314 #c2ff00\n0 315 #c4ff00\n0 316 #c7ff00\n0 317 #c9ff00\n0 318 #ccff00\n0 319 #ceff00\n0 320 #d1ff00\n0 321 #d3ff00\n0 322 #d6ff00\n0 323 #d8ff00\n0 324 #dbff00\n0 325 #ddff00\n0 326 #e0ff00\n0 327 #e2ff00\n0 328 #e5ff00\n0 329 #e7ff00\n0 330 #eaff00\n0 331 #ecff00\n0 332 #efff00\n0 333 #f1ff00\n0 334 #f4ff00\n0 335 #f6ff00\n0 336 #f9ff00\n0 337 #fbff00\n0 338 #feff00\n0 339 #fffe00\n0 340 #fffc00\n0 341 #fff900\n0 342 #fff700\n0 343 #fff400\n0 344 #fff200\n0 345 #ffef00\n0 346 #ffed00\n0 347 #ffea00\n0 348 #ffe800\n0 349 #ffe500\n0 350 #ffe300\n0 351 #ffe000\n0 352 #ffde00\n0 353 #ffdb00\n0 354 #ffd900\n0 355 #ffd600\n0 356 #ffd400\n0 357 #ffd100\n0 358 #ffcf00\n0 359 #ffcc00\n0 360 #ffca00\n0 361 #ffc700\n0 362 #ffc500\n0 363 #ffc200\n0 364 #ffc000\n0 365 #ffbd00\n0 366 #ffbb00\n0 367 #ffb800\n0 368 #ffb600\n0 369 #ffb300\n0 370 #ffb100\n0 371 #ffae00\n0 372 #ffac00\n0 373 #ffa900\n0 374 #ffa700\n0 375 #ffa400\n0 376 #ffa200\n0 377 #ff9f00\n0 378 #ff9d00\n0 379 #ff9a00\n0 380 #ff9800\n0 381 #ff9500\n0 382 #ff9300\n0 383 #ff9000\n0 384 #ff8e00\n0 385 #ff8b00\n0 386 #ff8900\n0 387 #ff8600\n0 388 #ff8400\n0 389 #ff8100\n0 390 #ff7f00\n0 391 #ff7c00\n0 392 #ff7a00\n0 393 #ff7700\n0 394 #ff7500\n0 395 #ff7200\n0 396 #ff7000\n0 397 #ff6d00\n0 398 #ff6b00\n0 399 #ff6800\n0 400 #ff6600\n0 401 #ff6300\n0 402 #ff6100\n0 403 #ff5e00\n0 404 #ff5c00\n0 405 #ff5900\n0 406 #ff5700\n0 407 #ff5400\n0 408 #ff5200\n0 409 #ff4f00\n0 410 #ff4d00\n0 411 #ff4a00\n0 412 #ff4800\n0 413 #ff4500\n0 414 #ff4300\n0 415 #ff4000\n0 416 #ff3e00\n0 417 #ff3b00\n0 418 #ff3900\n0 419 #ff3600\n0 420 #ff3400\n0 421 #ff3100\n0 422 #ff2f00\n0 423 #ff2c00\n0 424 #ff2a00\n0 425 #ff2700\n0 426 #ff2500\n0 427 #ff2200\n0 428 #ff2000\n0 429 #ff1d00\n0 430 #ff1b00\n0 431 #ff1800\n0 432 #ff1600\n0 433 #ff1300\n0 434 #ff1100\n0 435 #ff0e00\n0 436 #ff0c00\n0 437 #ff0900\n0 438 #ff0700\n0 439 #ff0400\n0 440 #ff0200\n0 441 #ff0000\n0 442 #ff0003\n0 443 #ff0005\n0 444 #ff0008\n0 445 #ff000a\n0 446 #ff000d\n0 447 #ff000f\n0 448 #ff0012\n0 449 #ff0014\n0 450 #ff0017\n0 451 #ff0019\n0 452 #ff001c\n0 453 #ff001e\n0 454 #ff0021\n0 455 #ff0023\n0 456 #ff0026\n0 457 #ff0028\n0 458 #ff002b\n0 459 #ff002d\n0 460 #ff0030\n0 461 #ff0032\n0 462 #ff0035\n0 463 #ff0037\n0 464 #ff003a\n0 465 #ff003c\n0 466 #ff003f\n0 467 #ff0041\n0 468 #ff0044\n0 469 #ff0046\n0 470 #ff0049\n0 471 #ff004b\n0 472 #ff004e\n0 473 #ff0050\n0 474 #ff0053\n0 475 #ff0055\n0 476 #ff0058\n0 477 #ff005a\n0 478 #ff005d\n0 479 #ff005f\n0 480 #ff0062\n0 481 #ff0064\n0 482 #ff0067\n0 483 #ff0069\n0 484 #ff006c\n0 485 #ff006e\n0 486 #ff0071\n0 487 #ff0073\n0 488 #ff0076\n0 489 #ff0078\n0 490 #ff007b\n0 491 #ff007d\n0 492 #ff0080\n0 493 #ff0082\n0 494 #ff0085\n0 495 #ff0087\n0 496 #ff008a\n0 497 #ff008c\n0 498 #ff008f\n0 499 #ff0091\n0 500 #ff0094\n0 501 #ff0096\n0 502 #ff0099\n0 503 #ff009b\n0 504 #ff009e\n0 505 #ff00a0\n0 506 #ff00a3\n0 507 #ff00a5\n0 508 #ff00a8\n0 509 #ff00aa\n0 510 #ff00ad\n0 511 #ff00af\n0 512 #ff00b2\n0 513 #ff00b4\n0 514 #ff00b7\n0 515 #ff00b9\n0 516 #ff00bc\n0 517 #ff00be\n0 518 #ff00c1\n0 519 #ff00c3\n0 520 #ff00c6\n0 521 #ff00c8\n0 522 #ff00cb\n0 523 #ff00cd\n0 524 #ff00d0\n0 525 #ff00d2\n0 526 #ff00d5\n0 527 #ff00d7\n0 528 #ff00da\n0 529 #ff00dc\n0 530 #ff00df\n0 531 #ff00e1\n0 532 #ff00e4\n0 533 #ff00e6\n0 534 #ff00e9\n0 535 #ff00eb\n0 536 #ff00ee\n0 537 #ff00f0\n0 538 #ff00f3\n0 539 #ff00f5\n0 540 #ff00f8\n0 541 #ff00fa\n0 542 #ff00fd\n0 543 #ff00ff\n";

/* Here are RNUM functions */
double randn();
long long int RGET(long long int *);
double R01GET(long long int *);
double RNGET(long long int *);
double RISIGET(long long int *,double);

/* Here are memory functions */
void meminit();
void memlink(void *,size_t);
void memdrop(void *);
void memprintf(int);
void memsearch(void *);
void * tmalloc(size_t);
void * tcalloc(size_t,size_t);
void * trealloc(void *,size_t);
void tfree(void *);

/* Here are llist functions */
struct litem * litemmake();
void litemtfree(struct litem * );
struct llist * llistmake();
struct llist * llistcopy();
struct llist * llistunion(struct llist *,struct llist *);
void llistgrowllist(struct llist *,struct llist *);
void llisttfree(struct llist *);
void llisttfree2(struct llist *);
void llisttfree3(struct llist *);
int isin(struct llist *,void *);
void litemadd(struct llist *,void *);
int litemexadd(struct llist *,void *);
void litemminus(struct llist *,void *);
struct litem * liteminsertafter(struct llist *,struct litem *);
struct litem * liteminsertbefore(struct llist *,struct litem *);
void llistkillfirst(struct llist *);
void llistkillfirst2(struct llist *);
void llistkilllast(struct llist *);
void llistkilllast2(struct llist *);
void llistprintf(struct llist *);
void llistprintf2(struct llist *);
int ra2ra_generic_compare_dictionary(void *,void *,void *);
int double_compare(void *,void *);
int int_compare(void *,void *);
int void_compare(void *,void *);
void llistsort(struct litem *,struct litem *,int,int (*)(void *,void *));
struct llist ** llistra2read(char *,int *);
int llistra2dump(struct llist **,int,char *);

/* Here are llitem functions */
struct llitem * llitemmake();
void llitemtfree(struct llitem *,void (*)(void *));
struct llitem * llitemaddorfind(int,struct llitem *,void *,int (*)(void *,void *));
struct llitem * llitemaddorfind_generic(int,struct llitem *,void *,int (*)(void *,void *,void *),void *);
struct llitem * llitemcoup(int,struct llitem *);
void llitembalance(struct llitem *);
struct llitem * llitemclimb(struct llitem *);
void llitemcheck(int,struct llitem *,int (*)(void *,void *));
void llitemprintf(struct llitem *,void (*)(void *));
int llitemlength(struct llitem *);
void llitemgrowllitem(struct llitem *,struct llitem *,int (*)(void *,void *));
void llistgrowllitem(struct llist *,struct llitem *);
void llist2llitem(struct llist *,struct llitem *,int (*)(void *,void *));
void llitem2llist(struct llitem *,struct llist *);

/* Here are the hist functions */
struct hist * histmake(int,double,double);
void histtfree(struct hist *);
void hist2file(void *,FILE *);
void * file2hist(FILE *,int *);
struct hist * histcopy(struct hist *);
void histadd(struct hist *,double,double);
void histprintf(struct hist *,char *);
void llist2hist(struct llist *,struct hist *);
void histdump(struct hist *,int,char *,char *,int);

/* Here are the statistics functions */
double znorm(double *z1r,double *z1i);
void zpz(double *,double *,double *,double *,double *,double *);
void zmz(double *,double *,double *,double *,double *,double *);
void zxz(double *,double *,double *,double *,double *,double *);
void zdz(double *,double *,double *,double *,double *,double *);
void zpeval(int,double *,double *,double *,double *,double *,double *);
int durand_kerner_rootfind(int,double *,double *,double *,double *);
int adi_round(double);
void doubleswap(double *,double *);
void periodify(char *,void *,void *,void *,void *);
void lliststats(struct llist *,double *,double *,double *,double *);
void lliststats2(struct llist *,double *,int,int);
double llistcorrelation(struct llist *,struct llist *);
void binstats(char *,void *,int,double *,double *);
int indexpack(int,int *,int *);
void indextract(int,int,int *,int *);
void maxmindex(char *,void *,int,int **,int *,int **,int *,double);
void stats(char *,void *,int,double *,double *,double *,double *);
void raaddra(double *,double *,int,int);
double radotmean(double *,double *,int);
double correlation(double *,double *,int);
double * spacesmear(double *,int,int,int);
void rgb2hsv(double,double,double,double *,double *,double *);
void hsv2rgb(double,double,double,double *,double *,double *);
void colorscale(int,double,double,double,double *,double *,double *);
void colorscaleinv(int,double,double,double,double,double,double *);
void fftwconvolve(void *,void *,double *,double *,double *,int,int);
double * ReadPNMfile(char *,int,int,int*,double *,double *,int*,int*);
int WritePNMfile(double *,int,int,double,double,char *);
int WritePNMfile_color(double *,int,int,double,double,char *,int);
void num2frame(int,char *);
int RescalePNMfiles(char *,int,int);
void rareset(void *,char *,int,void *);
void raprintf(void *,char *,int,int,char *);
void ra2jpg(void *,char *,int,int,int,char *,int);
void ra2jpg2(void *,char *,int,int,int,double,double,char *);
void * raread(char *,int *,int *,int *);
int checktofind_howmany(char *,int *,char *);
int checktofind(char *);
int multidradump(double *,int,int *,char *);
int multidraread(char *,int **,double **);
int radump(void *,char *,int,int,char *);
int rarecord(int,double *,int,char *);
int pnm2mpg(char *,int,int);
double ra_norm(double *,int);
double ra2ra_dot(double *,double *,int);
double * ra2ra_plus(double *,double *,int);
double * ra2ra_minus(double *,double *,int);
double * ra2ra_matrix_multiply(double *,int,int,int,double *,int,int,int);
void raplugin(double *,int,int,double *,int,int,int,int);
double * raplugout(double *,int,int,int,int,int,int);
void raplusequals(double *,int,double *);
void ratimesequals(double *,int,double);
void raaddequals(double *,int,double);
double * ra2power(void *,double *,int,double *,int,int);
void rara2corr(void *,void *,double **,int,int,double *,double *);

void ping();void pong();

