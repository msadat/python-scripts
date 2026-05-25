import streamlit as st
from datetime import date
from PIL import Image
import io
import base64

# -----------------------------
# Page Configuration
# -----------------------------

st.set_page_config(
    page_title="Executive Daily Construction Update",
    page_icon="🏗️",
    layout="wide"
)

# -----------------------------
# Helper Functions
# -----------------------------

def generate_executive_update(
    project_name,
    report_date,
    prepared_by,
    weather,
    overall_status,
    executive_summary,
    activities,
    quantities,
    issues,
    safety_notes,
    environmental_notes,
    upcoming_work,
    decisions_needed
):
    """Generate a polished executive-style daily construction update."""

    update_text = f"""
# {project_name} Daily Executive Update

**Date:** {report_date}  
**Prepared By:** {prepared_by}  
**Overall Status:** {overall_status}  

---

## Executive Summary

{executive_summary if executive_summary.strip() else "No executive summary provided."}

---

## Weather / Site Conditions

{weather if weather.strip() else "Weather/site conditions not provided."}

---

## Work Completed Today

{activities if activities.strip() else "No work activities entered."}

---

## Quantities / Production

{quantities if quantities.strip() else "No quantities entered."}

---

## Issues, Risks, or Delays

{issues if issues.strip() else "No major issues, risks, or delays reported."}

---

## Safety Notes

{safety_notes if safety_notes.strip() else "No safety incidents or concerns reported."}

---

## Environmental / SWPPP Notes

{environmental_notes if environmental_notes.strip() else "No environmental or SWPPP concerns reported."}

---

## Upcoming Work

{upcoming_work if upcoming_work.strip() else "No upcoming work entered."}

---

## Executive Decisions / Support Needed

{decisions_needed if decisions_needed.strip() else "No executive decisions or support required at this time."}

---

**End of Daily Update**
"""
    return update_text


def create_download_link(text, filename):
    """Create a downloadable text file link."""
    b64 = base64.b64encode(text.encode()).decode()
    href = f'<a href="data:file/txt;base64,{b64}" download="{filename}">Download Executive Update Text File</a>'
    return href


# -----------------------------
# App Header
# -----------------------------

st.title("🏗️ Executive Daily Construction Update Generator")
st.caption("Create professional daily construction updates with activities, quantities, risks, photos, and upcoming work.")

st.divider()

# -----------------------------
# Sidebar Inputs
# -----------------------------

st.sidebar.header("Report Information")

project_name = st.sidebar.text_input(
    "Project / Program Name",
    value="PA1490 SWEAT-1B"
)

report_date = st.sidebar.date_input(
    "Report Date",
    value=date.today()
)

prepared_by = st.sidebar.text_input(
    "Prepared By",
    value="Rafat Sadat"
)

overall_status = st.sidebar.selectbox(
    "Overall Status",
    [
        "On Track",
        "Minor Issues / Monitoring",
        "At Risk",
        "Delayed",
        "Critical Executive Attention Needed"
    ]
)

st.sidebar.divider()

st.sidebar.header("Report Format Options")

include_photos = st.sidebar.checkbox("Include Photo Section", value=True)
include_decision_section = st.sidebar.checkbox("Include Executive Decisions Section", value=True)

# -----------------------------
# Main Layout
# -----------------------------

tab1, tab2, tab3, tab4 = st.tabs(
    [
        "1. Daily Inputs",
        "2. Photos",
        "3. Generated Update",
        "4. Export"
    ]
)

# -----------------------------
# Tab 1: Daily Inputs
# -----------------------------

with tab1:
    st.header("Daily Construction Inputs")

    col1, col2 = st.columns(2)

    with col1:
        weather = st.text_area(
            "Weather / Site Conditions",
            height=120,
            placeholder="Example: Clear weather, normal site access, no weather-related impacts."
        )

        executive_summary = st.text_area(
            "Executive Summary",
            height=160,
            placeholder=(
                "Example: Work progressed as planned today with continued excavation on TWY C inner ARFF road, "
                "ongoing ASF building installation, and fence work on the west airfield. No major safety issues were reported."
            )
        )

        activities = st.text_area(
            "Work Completed Today",
            height=250,
            placeholder=(
                "Example:\n"
                "- TWY C inner ARFF road excavation continued.\n"
                "- Removed old bridge outside AOA at old W. ARFF road.\n"
                "- Installed geotextile for outer ARFF road.\n"
                "- ASF building crane operation completed; installed 3 arches and 2 flat ends."
            )
        )

    with col2:
        quantities = st.text_area(
            "Quantities / Production",
            height=180,
            placeholder=(
                "Example:\n"
                "- 150 LF excavation for triple barrel drainage structure.\n"
                "- 150 ft x 60 ft subgrade preparation completed.\n"
                "- 600 CY unclassified excavation hauled off.\n"
                "- 500 LF temporary AOA fence removed."
            )
        )

        issues = st.text_area(
            "Issues, Risks, or Delays",
            height=160,
            placeholder=(
                "Example:\n"
                "- No major delays reported.\n"
                "- Continue monitoring creek work sequencing pending aquatic relocation and nest survey."
            )
        )

        safety_notes = st.text_area(
            "Safety Notes",
            height=120,
            placeholder="Example: No safety incidents reported. Crews maintained AOA access control requirements."
        )

        environmental_notes = st.text_area(
            "Environmental / SWPPP Notes",
            height=120,
            placeholder="Example: No SWPPP deficiencies reported. Creek work remains restricted pending required environmental clearance."
        )

        upcoming_work = st.text_area(
            "Upcoming Work",
            height=160,
            placeholder=(
                "Example:\n"
                "- Concrete crew to continue forming TWY E.\n"
                "- Scheduled 180 CY hand pour tomorrow.\n"
                "- Continue storm pipe construction and ARFF road preparation."
            )
        )

        if include_decision_section:
            decisions_needed = st.text_area(
                "Executive Decisions / Support Needed",
                height=120,
                placeholder="Example: No executive decisions required at this time."
            )
        else:
            decisions_needed = ""

# -----------------------------
# Tab 2: Photos
# -----------------------------

uploaded_photos = []
photo_captions = []

with tab2:
    st.header("Daily Photos")

    if include_photos:
        st.info("Upload daily progress photos and add a short caption for each photo.")

        photo_files = st.file_uploader(
            "Upload Progress Photos",
            type=["jpg", "jpeg", "png"],
            accept_multiple_files=True
        )

        if photo_files:
            for idx, photo_file in enumerate(photo_files):
                st.subheader(f"Photo {idx + 1}")

                image = Image.open(photo_file)
                st.image(image, use_container_width=True)

                caption = st.text_input(
                    f"Caption for Photo {idx + 1}",
                    key=f"caption_{idx}",
                    placeholder="Example: TWY C inner ARFF road excavation progress."
                )

                uploaded_photos.append(photo_file)
                photo_captions.append(caption)

                st.divider()
        else:
            st.warning("No photos uploaded yet.")
    else:
        st.warning("Photo section is turned off in the sidebar.")

# -----------------------------
# Generate Report Text
# -----------------------------

generated_update = generate_executive_update(
    project_name=project_name,
    report_date=report_date.strftime("%B %d, %Y"),
    prepared_by=prepared_by,
    weather=weather,
    overall_status=overall_status,
    executive_summary=executive_summary,
    activities=activities,
    quantities=quantities,
    issues=issues,
    safety_notes=safety_notes,
    environmental_notes=environmental_notes,
    upcoming_work=upcoming_work,
    decisions_needed=decisions_needed
)

# Add photo captions to generated update
if include_photos:
    generated_update += "\n\n---\n\n## Progress Photos\n\n"

    if photo_captions:
        for idx, caption in enumerate(photo_captions):
            generated_update += f"**Photo {idx + 1}:** {caption if caption.strip() else 'No caption provided.'}\n\n"
    else:
        generated_update += "No photos uploaded.\n"

# -----------------------------
# Tab 3: Generated Update
# -----------------------------

with tab3:
    st.header("Generated Executive Update")

    st.markdown(generated_update)

    st.divider()

    st.subheader("Copy-Friendly Version")
    st.text_area(
        "Copy this text into email, Teams, or a daily report.",
        value=generated_update,
        height=500
    )

# -----------------------------
# Tab 4: Export
# -----------------------------

with tab4:
    st.header("Export Daily Update")

    filename = f"{project_name.replace(' ', '_')}_Daily_Update_{report_date.strftime('%Y_%m_%d')}.txt"

    st.markdown(
        create_download_link(generated_update, filename),
        unsafe_allow_html=True
    )

    st.divider()

    st.subheader("Suggested Email Subject")
    subject_line = f"{project_name} Daily Executive Update - {report_date.strftime('%m/%d/%Y')}"
    st.code(subject_line)

    st.subheader("Suggested Executive Email Opening")
    email_opening = f"""
Good evening,

Please see below for the {project_name} daily executive update for {report_date.strftime('%B %d, %Y')}.

Overall status: {overall_status}.

"""
    st.text_area(
        "Email Opening",
        value=email_opening,
        height=140
    )

# -----------------------------
# Footer
# -----------------------------

st.divider()
st.caption("Daily update generator developed for construction management, airfield programs, and executive reporting.")
