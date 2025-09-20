#include "tables.h"

#include "types.h"
#include <iostream>

#include <cstring>

//Reverses a bitboard                        
Bitboard reverse(Bitboard b) {
	b = (b & 0x5555555555555555) << 1 | (b >> 1) & 0x5555555555555555;
	b = (b & 0x3333333333333333) << 2 | (b >> 2) & 0x3333333333333333;
	b = (b & 0x0f0f0f0f0f0f0f0f) << 4 | (b >> 4) & 0x0f0f0f0f0f0f0f0f;
	b = (b & 0x00ff00ff00ff00ff) << 8 | (b >> 8) & 0x00ff00ff00ff00ff;

	return (b << 48) | ((b & 0xffff0000) << 16) |
		((b >> 16) & 0xffff0000) | (b >> 48);
}

//Calculates sliding attacks from a given square, on a given axis, taking into
//account the blocking pieces. This uses the Hyperbola Quintessence Algorithm.
Bitboard sliding_attacks(Square square, Bitboard occ, Bitboard mask) {
	return (((mask & occ) - SQUARE_BB[square] * 2) ^
		reverse(reverse(mask & occ) - reverse(SQUARE_BB[square]) * 2)) & mask;
}

//Returns rook attacks from a given square, using the Hyperbola Quintessence Algorithm. Only used to initialize
//the magic lookup table
Bitboard get_rook_attacks_for_init(Square square, Bitboard occ) {
	return sliding_attacks(square, occ, MASK_FILE[file_of(square)]) |
		sliding_attacks(square, occ, MASK_RANK[rank_of(square)]);
}

//Initializes the magic lookup table for rooks
void initialise_rook_attacks() {
	Bitboard edges, subset, index;

	for (Square sq = a1; sq <= h8; ++sq) {
		edges = ((MASK_RANK[AFILE] | MASK_RANK[HFILE]) & ~MASK_RANK[rank_of(sq)]) |
			((MASK_FILE[AFILE] | MASK_FILE[HFILE]) & ~MASK_FILE[file_of(sq)]);
		ROOK_ATTACK_MASKS[sq] = (MASK_RANK[rank_of(sq)]
			^ MASK_FILE[file_of(sq)]) & ~edges;
		ROOK_ATTACK_SHIFTS[sq] = 64 - pop_count(ROOK_ATTACK_MASKS[sq]);

		subset = 0;
		do {
			index = subset;
			index = index * ROOK_MAGICS[sq];
			index = index >> ROOK_ATTACK_SHIFTS[sq];
			ROOK_ATTACKS[sq][index] = get_rook_attacks_for_init(sq, subset);
			subset = (subset - ROOK_ATTACK_MASKS[sq]) & ROOK_ATTACK_MASKS[sq];
		} while (subset);
	}
}

//Returns the attacks bitboard for a rook at a given square, using the magic lookup table
constexpr Bitboard get_rook_attacks(Square square, Bitboard occ) {
	return ROOK_ATTACKS[square][((occ & ROOK_ATTACK_MASKS[square]) * ROOK_MAGICS[square])
		>> ROOK_ATTACK_SHIFTS[square]];
}

//Returns the 'x-ray attacks' for a rook at a given square. X-ray attacks cover squares that are not immediately
//accessible by the rook, but become available when the immediate blockers are removed from the board 
Bitboard get_xray_rook_attacks(Square square, Bitboard occ, Bitboard blockers) {
	Bitboard attacks = get_rook_attacks(square, occ);
	blockers &= attacks;
	return attacks ^ get_rook_attacks(square, occ ^ blockers);
}

//Returns bishop attacks from a given square, using the Hyperbola Quintessence Algorithm. Only used to initialize
//the magic lookup table
Bitboard get_bishop_attacks_for_init(Square square, Bitboard occ) {
	return sliding_attacks(square, occ, MASK_DIAGONAL[diagonal_of(square)]) |
		sliding_attacks(square, occ, MASK_ANTI_DIAGONAL[anti_diagonal_of(square)]);
}

//Initializes the magic lookup table for bishops
void initialise_bishop_attacks() {
	Bitboard edges, subset, index;

	for (Square sq = a1; sq <= h8; ++sq) {
		edges = ((MASK_RANK[AFILE] | MASK_RANK[HFILE]) & ~MASK_RANK[rank_of(sq)]) |
			((MASK_FILE[AFILE] | MASK_FILE[HFILE]) & ~MASK_FILE[file_of(sq)]);
		BISHOP_ATTACK_MASKS[sq] = (MASK_DIAGONAL[diagonal_of(sq)]
			^ MASK_ANTI_DIAGONAL[anti_diagonal_of(sq)]) & ~edges;
		BISHOP_ATTACK_SHIFTS[sq] = 64 - pop_count(BISHOP_ATTACK_MASKS[sq]);

		subset = 0;
		do {
			index = subset;
			index = index * BISHOP_MAGICS[sq];
			index = index >> BISHOP_ATTACK_SHIFTS[sq];
			BISHOP_ATTACKS[sq][index] = get_bishop_attacks_for_init(sq, subset);
			subset = (subset - BISHOP_ATTACK_MASKS[sq]) & BISHOP_ATTACK_MASKS[sq];
		} while (subset);
	}
}

//Returns the attacks bitboard for a bishop at a given square, using the magic lookup table
constexpr Bitboard get_bishop_attacks(Square square, Bitboard occ) {
	return BISHOP_ATTACKS[square][((occ & BISHOP_ATTACK_MASKS[square]) * BISHOP_MAGICS[square])
		>> BISHOP_ATTACK_SHIFTS[square]];
}

//Returns the 'x-ray attacks' for a bishop at a given square. X-ray attacks cover squares that are not immediately
//accessible by the rook, but become available when the immediate blockers are removed from the board 
Bitboard get_xray_bishop_attacks(Square square, Bitboard occ, Bitboard blockers) {
	Bitboard attacks = get_bishop_attacks(square, occ);
	blockers &= attacks;
	return attacks ^ get_bishop_attacks(square, occ ^ blockers);
}

//Initializes the lookup table for the bitboard of squares in between two given squares (0 if the 
//two squares are not aligned)
void initialise_squares_between() {
	Bitboard sqs;
	for (Square sq1 = a1; sq1 <= h8; ++sq1)
		for (Square sq2 = a1; sq2 <= h8; ++sq2) {
			sqs = SQUARE_BB[sq1] | SQUARE_BB[sq2];
			if (file_of(sq1) == file_of(sq2) || rank_of(sq1) == rank_of(sq2))
				SQUARES_BETWEEN_BB[sq1][sq2] =
				get_rook_attacks_for_init(sq1, sqs) & get_rook_attacks_for_init(sq2, sqs);
			else if (diagonal_of(sq1) == diagonal_of(sq2) || anti_diagonal_of(sq1) == anti_diagonal_of(sq2))
				SQUARES_BETWEEN_BB[sq1][sq2] =
				get_bishop_attacks_for_init(sq1, sqs) & get_bishop_attacks_for_init(sq2, sqs);
		}
}


//Initializes the lookup table for the bitboard of all squares along the line of two given squares (0 if the 
//two squares are not aligned)
void initialise_line() {
	for (Square sq1 = a1; sq1 <= h8; ++sq1)
		for (Square sq2 = a1; sq2 <= h8; ++sq2) {
			if (file_of(sq1) == file_of(sq2) || rank_of(sq1) == rank_of(sq2))
				LINE[sq1][sq2] =
				get_rook_attacks_for_init(sq1, 0) & get_rook_attacks_for_init(sq2, 0)
				| SQUARE_BB[sq1] | SQUARE_BB[sq2];
			else if (diagonal_of(sq1) == diagonal_of(sq2) || anti_diagonal_of(sq1) == anti_diagonal_of(sq2))
				LINE[sq1][sq2] =
				get_bishop_attacks_for_init(sq1, 0) & get_bishop_attacks_for_init(sq2, 0)
				| SQUARE_BB[sq1] | SQUARE_BB[sq2];
		}
}

//Initializes the table containg pseudolegal attacks of each piece for each square. This does not include blockers
//for sliding pieces
void initialise_pseudo_legal() {
	memcpy(PAWN_ATTACKS[WHITE], WHITE_PAWN_ATTACKS, sizeof(WHITE_PAWN_ATTACKS));
	memcpy(PAWN_ATTACKS[BLACK], BLACK_PAWN_ATTACKS, sizeof(BLACK_PAWN_ATTACKS));
	memcpy(PSEUDO_LEGAL_ATTACKS[KNIGHT], KNIGHT_ATTACKS, sizeof(KNIGHT_ATTACKS));
	memcpy(PSEUDO_LEGAL_ATTACKS[KING], KING_ATTACKS, sizeof(KING_ATTACKS));
	for (Square s = a1; s <= h8; ++s) {
		PSEUDO_LEGAL_ATTACKS[ROOK][s] = get_rook_attacks_for_init(s, 0);
		PSEUDO_LEGAL_ATTACKS[BISHOP][s] = get_bishop_attacks_for_init(s, 0);
		PSEUDO_LEGAL_ATTACKS[QUEEN][s] = PSEUDO_LEGAL_ATTACKS[ROOK][s] |
			PSEUDO_LEGAL_ATTACKS[BISHOP][s];
	}
}

//Initializes lookup tables for rook moves, bishop moves, in-between squares, aligned squares and pseudolegal moves
void initialise_all_databases() {
	initialise_rook_attacks();
	initialise_bishop_attacks();
	initialise_squares_between();
	initialise_line();
	initialise_pseudo_legal();
}
