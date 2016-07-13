/*
 * ctfutils.h - declarations for utility functions used by ctfplotter and
 *                 ctfphaseflip
 *
 *  $Id: ctfutils.h,v d8d1ce772290 2010/04/02 00:03:31 mast $
 */
#ifndef CTFUTILS_H
#define CTFUTILS_H

typedef struct ilist_struct Ilist;

typedef struct defocus_struct {
  int startingSlice;
  int endingSlice;
  double lAngle;
  double hAngle;
  double defocus;
} SavedDefocus;

float *readTiltAngles(const char *angleFile, int mNzz, float mAngleSign,
                      float &mMinAngle, float &mMaxAngle);
void addItemToDefocusList(Ilist *mSaved, SavedDefocus toSave);
Ilist *readDefocusFile(const char *mFnDefocus);
int checkAndFixDefocusList(Ilist *list, float *angles, int nz);

#endif
