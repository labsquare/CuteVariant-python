from PySide2.QtWidgets import QApplication

__title__ = "Columns selector"
__description__ = "A plugin to display fields "
__long_description__ = QApplication.instance().translate(
    "fields_editor",
    """
<p>This plugin is used to manage the fields visible from the main "variant view"
window.</p>

The fields are divided into 3 main categories:
<ul>
<li><b>variants:</b> fields concerning the variants in their globality
(chromosome, position, reference base, alternative base, etc.).</li>
<li><b>annotations:</b> fields added by annotation tools such as snpeff (gene,
consequence of the variant, impact rating, etc.).</li>
<li><b>samples:</b> fields concerning the individuals sequenced in the project.</li>
</ul>

<p>As the number of fields is large, it may be interesting to <em>use the keyword
search button</em> to quickly find those that are relevant in a study.</p>
""",
)
__author__ = "Sacha schutz"
__version__ = "1.0.0"
